#include "mbconfig.h"
#include "userel.h"
#include "elem.h"
#include "dataman.h"
#include "dofown.h"
#include "strnode.h"
#include "matvec.h"
#include "varsubmat.h"
#include "varsubvec.h"
#include <vector>

static const doublereal gVal = 9.80665;
static const doublereal epsFD = 1e-6;

// Elem と UserDefinedElem の二つの基底クラスから継承
// virtual → UserDefinedElem も内部で Elem を継承している
class ModuleLM : virtual public Elem, public UserDefinedElem {
    
    unsigned Seg; // セグメント数
    doublereal EA, CA; // EA:軸方向剛性, CA:軸方向減衰係数
    doublereal zBed, kn, cn, mu; // 海底の垂直位置，海底の法線方向剛性，海底の法線方向減衰，摩擦係数
    
    // std::vector 動的配列
    std::vector<doublereal> L0; // 各セグメントの初期長さ
    std::vector<doublereal> M; // 各セグメントの質量
    
    // StructNode * ポインタ型でノードオブジェクトへの参照
    std::vector<StructNode *> N; // 構造ノードへのポインタ配列

public:
    
    // uL:要素のラベル番号，DofOwner *pDO: 自由度の所有者オブジェクト（読み取り専用），DataManager *pDM: データ管理オブジェクト（更新可能），MBDynParser &HP: パーサーオブジェクト（参照渡し）
    // Elem と UserDefinedElem の両方を初期化
    ModuleLM(unsigned uL, const DofOwner *pDO, DataManager *pDM, MBDynParser &HP) : Elem(uL, pDO), UserDefinedElem(uL) {
        
        unsigned fair = HP.GetInt(); // FP のノード番号の取得
        unsigned anchor = HP.GetInt(); // AP のノード番号の取得
        
        // パラメータ
        if (HP.IsKeyWord("seg")) Seg = HP.GetInt(); else Seg = 20;
        if (HP.IsKeyWord("EA")) EA = HP.GetReal(); else EA = 1.e8;
        if (HP.IsKeyWord("CA")) CA = HP.GetReal(); else CA = 5.e5;
        doublereal l0 = 45.; // 初期長さ
        doublereal m0 = 120.; // 質量
        if (HP.IsKeyWord("l0")) l0 = HP.GetReal();
        if (HP.IsKeyWord("mass")) m0 = HP.GetReal();
        
        // assign → 全要素を同じ値で埋める
        L0.assign(Seg, l0); // 全セグメントの長さを同じに設定
        M.assign(Seg - 1, m0); // 内部ノードの質量を設定
        
        // 海底境界条件
        zBed = -2.e2; // 海底位置 -200m
        kn = 1.e7; // 法線剛性
        cn = 1.e5; // 法線減衰
        mu = 0.6; // 摩擦係数
        
        if (HP.IsKeyWord("seabed")) {
            zBed = HP.GetReal(); kn = HP.GetReal(); cn = HP.GetReal(); mu = HP.GetReal();
        }

        // Seg + 1個の要素を確保（両端＋中間ノード）→ 要素Seg個ならノードはSeg + 1個必要
        N.resize(Seg + 1); // 係留索の各分割点用のノード配列を準備
        N[0] = dynamic_cast<StructNode *>(pDM->ReadNode(fair)); // 境界ノードの設定：浮体側
        N[Seg] = dynamic_cast<StructNode *>(pDM->ReadNode(anchor)); // 境界ノードの設定：係留店側
        for (unsigned i = 1; i < Seg; ++i) {
            Vec3 X0 = N[0]->GetXCurr() + (N[Seg]->GetXCurr() - N[0]->GetXCurr()) * (static_cast<doublereal>(i) / Seg); // 両端のノード間を線形補間：i/Seg で 0 から 1 の線形分布を生成
            StructNode *pN = pDM->CreateStructNode(pDO, X0, Vec3(0., 0., 0.)); // 初期速度 0 の新しいノードの作成
            N[i] = pN;
        }
    }

    // 基底クラスへのポインタを通じた動的削除時に適切なデストラクタが呼ばれる
    virtual ~ModuleLM() {}
    // この要素が Mooring Line であることを明示
    virtual std::string GetElemType() const { return "ModuleMooringLM"; }
    // 両端は境界条件として固定
    // 自由度の管理：内部ノード数：Seg - 1（両端を除く）, 各ノード：3自由度（x, y, z座標）, 合計：(Seg - 1) × 3 自由度 → 座標のみで回転自由度はなし
    // 両端のーどは既存要素に固定されているが，内部ノードは未知数
    virtual unsigned int iGetNumDof() const { return (Seg - 1) * 3; }
    // 行列サイズの設定 → 残差ベクトルのサイズなどに重要
    // メモリの解放
    virtual void WorkSpaceDim(integer *piNumRows, integer *piNumCols) const {
        *piNumRows = *piNumCols = (Seg - 1) * 3; // 行列次元は自由度数と一致
    }

    // const → 読み取り専用の関数
    Vec3 Tension(unsigned i, unsigned j) const {
        
        // rij: ノード i から j への位置ベクトル
        // ->: ポインタを通したメンバーアクセス
        Vec3 rij = N[j]->GetXCurr() - N[i]->GetXCurr();
        // セグメントの現在長
        doublereal dist = rij.Norm();
        if (dist == 0.) return Vec3(0., 0., 0.); // 距離が 0 の場合の回避
        // 単位方向ベクトル
        Vec3 t = rij / dist;
        
        // 材料の力学特性
        doublereal dL = dist - L0[i]; // dL: セグメントの伸び量, L0[i]: 初期長
        doublereal k = EA / L0[i]; // 軸方向剛性を長さで正規化
        doublereal c = CA / L0[i]; // 軸方向減衰係数を長さで正規化
        doublereal vrel = Dot(N[j]->GetVCurr() - N[i]->GetVCurr(), t); // 2ノード間の相対速度の軸方向成分

        // 張力の計算
        // 弾性力 k * dL, 減衰力 c * vrel, 単位ベクトル t
        return t * (k * dL + c * vrel);
    }

    Vec3 SeabedForce(unsigned i) const {
        
        const Vec3 &x = N[i]->GetXCurr();
        const Vec3 &v = N[i]->GetVCurr();
        
        // z 座標を取得し，海底より上なら接触力なし
        if (x(3) >= zBed) return Vec3(0., 0., 0.);

        // 垂直方向の接触力（法線力）
        doublereal dn = zBed - x(3); // 海底への貫入量
        doublereal fn = kn * dn - cn * v(3); // 弾性力: kn * dn（剛性 × 貫入深さ）, 減衰力: cn * v(3)（減衰係数 × 垂直速度）
        Vec3 fnVec(0., 0., fn); // 垂直方向のみの力

        // 水平方向の摩擦力
        Vec3 vt = v - Vec3(0., 0., v(3)); // 垂直成分を除去して vt: 水平方向の速度ベクトルを取得
        doublereal vtMag = vt.Norm();
        if (vtMag > 0.) {
            // クーロン摩擦
            // 方向: 速度と逆方向（-vt / vtMag）
            // 大きさ: 静止摩擦力と動摩擦力の最小値にする？
            Vec3 ft = -vt / vtMag * mu * std::abs(fn); // 静的クーロン摩擦
            return fnVec + ft;
        }
        return fnVec;
    }

    // AssRes 関数により係留索の各内部ノードにおける動的つり合い式が計算され，MBDyn の時間積分スキームに提供される
    virtual SubVectorHandler &AssRes(SubVectorHandler &RH, doublereal dCoef, const VectorHandler &, const VectorHandler &) {
        
        // 内部ノードの自由度数に設定し，0 で初期化
        RH.SetSizeZero((Seg - 1) * 3);
        
        for (unsigned i = 1; i < Seg; ++i) {
            // 運動方程式
            // 前のセグメントからの張力, 次のセグメントからの張力, 重力, 海底反力と摩擦力, 慣性力(質量 × 加速度)
            // R = F_external - M * a が 0 になる時，動的平衡状態にある
            Vec3 F = Tension(i - 1, i) + Tension(i + 1, i) + Vec3(0., 0., -M[i - 1] * gVal) + SeabedForce(i) - M[i - 1] * N[i]->GetXPrimePrime();
            // 残差の格納
            RH.Put(3 * (i - 1) + 1, F(1)); // x 成分
            RH.Put(3 * (i - 1) + 2, F(2)); // y 成分
            RH.Put(3 * (i - 1) + 3, F(3)); // z 成分
        }
        return RH;
    }

    // WH: 行列計算結果の格納先, XCurr: 現在の位置ベクトル, XPrimeCurr: 現在の速度ベクトル
    // AssRes 関数で計算した残差ベクトルを状態ベクトルに関して微分したもの
    virtual VariableSubMatrixHandler &AssJac(VariableSubMatrixHandler &WH, const VectorHandler &XCurr, const VectorHandler &XPrimeCurr) {
        
        // 密行列として設定し，全要素を 0 で初期化
        const integer dim = (Seg - 1) * 3;
        DenseSubMatrixHandler &K = WH.SetDense(dim, dim);
        K.Zero();

        // 基準残差の計算
        // 現在状態での残差ベクトル
        ColumnVector R0(dim);
        SubVectorHandler RH(dim);
        AssRes(RH, 0., XCurr, XPrimeCurr);

        // 配列のコピー
        for (integer i = 1; i <= dim; ++i) R0(i) = RH(i);

        // 数値微分
        for (integer k = 1; k <= dim; ++k) {
            // XCurr の複製を作成
            VectorHandler Xpert = XCurr;
            // 一つの成分のみ変更，各自由度を小さく摂動
            // ある自由度（座標の 3 自由度のどれか）をほんの少し動かして，システム全体の応答（残差）がどう変わるかを見る
            Xpert.Put(k, Xpert(k) + epsFD);
            // 摂動を加えた R1 ベクトル
            SubVectorHandler R1(dim);
            // 摂動後の残差ベクトルの計算
            AssRes(R1, 0., Xpert, XPrimeCurr);
            // 差分で偏微分を近似
            for (integer r = 1; r <= dim; ++r) K(r, k) = (R1(r) - R0(r)) / epsFD;
            // 数値安定化処理
            K(k, k) += 1.e3;
        }
        return WH;
    }

    // 出力
    virtual void Output(OutputHandler &OH) const {
        // 要素ラベル，ノード番号，位置ベクトルを出力
        for (unsigned i = 0; i <= Seg; ++i) OH << GetLabel() << " " << i << " " << N[i]->GetXCurr() << std::endl;
    }

    // リスタート機能：何も実装していない
    virtual std::ostream &Restart(std::ostream &out) const { return out; }
    // プライベートデータ出力：何も実装していない
    virtual doublereal dGetPrivData(unsigned int) const { return 0.; }
};

// 動的要素の生成
static UserDefinedElem *CreateModuleLM(unsigned uLabel, const DofOwner *pDO, DataManager *pDM, MBDynParser &HP) {
    return new ModuleLM(uLabel, pDO, pDM, HP);
}

// MBDyn への要素登録
static UserDefElemCreator regModuleLM("ModuleMooringLM", CreateModuleLM);
