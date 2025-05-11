/*  module-lm.cc ─ Lumped-mass mooring element with seabed reaction (MBDyn 1.7.3)  */
/*  ────────────────────────────────────────────────────────────────────────────── */
#include "mbconfig.h"
#include "elem.h"
#include "node.h"
#include "dataman.h"
#include "userdefinedelem.h"
#include "solver.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>
#include <memory>

/* === 任意海底形状を扱う基底クラス ========================================== */
class Seabed {
public:
    virtual doublereal GetZ(doublereal x, doublereal y) const = 0;
    virtual ~Seabed() {}
};
/* 定数深度 (flat) */
class FlatSeabed final : public Seabed {
    const doublereal z0;
public:
    explicit FlatSeabed(doublereal z) : z0(z) {}
    doublereal GetZ(doublereal, doublereal) const override { return z0; }
};
/* 格子ファイル (x y z) ASCII */
class GridSeabed final : public Seabed {
    struct Cell { doublereal x,y,z; };
    std::vector<Cell> grid;
public:
    explicit GridSeabed(const std::string& fname)
    {
        std::ifstream fin(fname.c_str());
        if(!fin) silent_cerr("GridSeabed: cannot open "<<fname<<std::endl);
        Cell c;
        while(fin>>c.x>>c.y>>c.z) grid.push_back(c);
    }
    doublereal GetZ(doublereal x, doublereal y) const override
    {
        doublereal dz = std::numeric_limits<doublereal>::max();
        doublereal z  = 0.;
        for(const auto& c: grid){
            doublereal d2 = sqr(c.x-x)+sqr(c.y-y);
            if(d2<dz){ dz=d2; z=c.z; }
        }
        return z;                         /* 最近点の z を返す簡易実装          */
    }
};

/* === 係留要素 =============================================================== */
class ModuleLM final : virtual public Elem,
                       public UserDefinedElem {
    /* ─ 入力パラメータ ─ */
    const unsigned Seg;
    const doublereal g;
    const doublereal kn, cn, mu;          /* seabed parameters                */

    /* 区間毎性質（将来: 中央だけ繊維など） */
    struct Prop { doublereal l0,m,EA; };
    std::vector<Prop> P;                  /* size = Seg                       */

    /* ノード配列 */
    std::vector<Node*> N;                 /* N[0]…N[Seg]                      */

    /* 海底プロファイル */
    std::unique_ptr<Seabed> pBed;

public:
    /* ───────────────────────── コンストラクタ ───────────────────────── */
    ModuleLM(unsigned uLabel,
             const DofOwner* pDO,
             DataManager* pDM,
             MBDynParser& HP)
    : Elem(uLabel,pDO),
      UserDefinedElem(uLabel),
      Seg(HP.GetInt("segments")),
      g  (HP.GetRealOrDefault("gravity",9.80665)),
      kn (HP.GetRealOrDefault("seabed kn",1.e7)),
      cn (HP.GetRealOrDefault("seabed cn",1.e5)),
      mu (HP.GetRealOrDefault("seabed mu",0.6))
    {
        /* ─ 区間物性読み取り ─ */
        P.resize(Seg);
        if(HP.IsKeyWord("property array")){
            for(unsigned i=0;i<Seg;i++){
                HP.GetLine();                        /* 格納行: l0 m EA */
                P[i].l0 = HP.GetReal();
                P[i].m  = HP.GetReal();
                P[i].EA = HP.GetReal();
            }
        }else{
            doublereal l0 = HP.GetReal("l0");
            doublereal m  = HP.GetReal("mass");
            doublereal EA = HP.GetReal("EA");
            for(auto& p:P){ p = {l0,m,EA}; }
        }

        /* ─ ノード ─ */
        for(unsigned i=0;i<=Seg;i++){
            Node* pN = 0;
            if(HP.IsKeyWord("existing_node")){
                pN = pDM->FindNode(HP.GetInt("existing_node"));
            }else{
                pN = pDM->ReadNode(Node::STRUCTURAL);
            }
            N.push_back(pN);
            AddNode(pN);
        }

        /* ─ 海底プロファイル ─ */
        if(HP.IsKeyWord("seabed file")){
            pBed.reset(new GridSeabed(HP.GetString()));
        }else{
            doublereal z0 = HP.GetRealOrDefault("seabed z0",0.);
            pBed.reset(new FlatSeabed(z0));
        }
    }

    ~ModuleLM() override = default;

    /* ========== MBDyn インタフェース最小実装 =============================== */
    unsigned    iGetNumDof()           const override { return 0; }
    unsigned    iGetNumPrivData()      const override { return 0; }
    doublereal  dGetPrivData(unsigned) const override { return 0.; }
    unsigned    iGetNumConnectedNodes()const override { return Seg+1; }
    const Node* pGetNode(unsigned i)   const override { return N[i]; }
    unsigned    iGetInitialNumDof()    const override { return 0; }
    std::ostream& Restart(std::ostream& o) const override { return o; }

    /* --- 初期形状 -------------------------------------------------------- */
    void SetInitialValue(VectorHandler& X, VectorHandler&) override
    {
        const Vec3 p0 = N[0]->GetXCurr();
        const Vec3 pL = N[Seg]->GetXCurr();
        const Vec3 d  = (pL-p0)/Seg;
        for(unsigned i=1;i<Seg;i++){
            N[i]->PutXCurr(p0+doublereal(i)*d);
            N[i]->PutVCurr(Zero3);
        }
    }

    /* --- Future stub functions ------------------------------------------ */
    void ApplyHydrodynamicForce(unsigned i, Vec3& Fi) const
    {
        /* TODO: 流体力（静水圧・Morison 等）を計算して Fi へ加算 */
    }
    void UpdateWaveKinematics(doublereal t)
    {
        /* TODO: 波浪場による流速・加速度を内部状態に更新 */
    }
    Prop GetMaterial(unsigned seg) const
    {
        /* TODO: seg が合成繊維区間なら別素材を返す  */
        return P[seg];
    }

    /* --- 残差 ----------------------------------------------------------- */
    SubVectorHandler&
    AssRes(SubVectorHandler& RH, doublereal dCoef,
           const VectorHandler&, const VectorHandler&) override
    {
        UpdateWaveKinematics(pGetDofOwner()->dGetTime());     /* 将来用 */
        const unsigned int dim=(Seg-1)*3+2*3;
        RH.SetSize(GetFirstIndex(),GetFirstIndex()+dim-1); RH.Zero();

        for(unsigned i=0;i<=Seg;i++){
            Vec3 Fi = Zero3;

            /* ─ 張力 ─ */
            if(i>0){
                const Prop left = GetMaterial(i-1);
                Vec3 dL = N[i]->GetXCurr()-N[i-1]->GetXCurr();
                doublereal LL = dL.Norm();
                Vec3 TL = left.EA*(LL-left.l0)/LL * dL;
                Fi += TL;
                if(i<Seg) RH.Put(GetFirstIndex()+(i-1)*3+0,-TL(1));
            }
            if(i<Seg){
                const Prop right = GetMaterial(i);
                Vec3 dR = N[i+1]->GetXCurr()-N[i]->GetXCurr();
                doublereal LR = dR.Norm();
                Vec3 TR = right.EA*(LR-right.l0)/LR * dR;
                Fi -= TR;
            }

            /* ─ 慣性 ─ */
            if(i>0 && i<Seg){
                Fi += GetMaterial(i).m * dCoef * N[i]->GetVCurr();
                RH.Put(GetFirstIndex()+(i-1)*3+0,-Fi(1));
            }

            /* ─ 海底反力・摩擦 ─ */
            const Vec3 pos = N[i]->GetXCurr();
            const Vec3 vel = N[i]->GetVCurr();
            const doublereal zb = pBed->GetZ(pos(1),pos(2));          /* 底面 z */
            if(pos(3)<zb){                                            /* 浸入時 */
                const doublereal pen = zb-pos(3);
                const doublereal vrel = -vel(3);                      /* 上向き正 */
                const doublereal Fn = kn*pen + cn*vrel;
                Vec3 FnVec(0.,0.,Fn);
                Fi += FnVec;

                /* Coulomb friction */
                Vec3 vxy(vel(1),vel(2),0.);
                doublereal vabs = vxy.Norm();
                if(vabs>1e-6){
                    Vec3 Ft = -mu*Fn * vxy/vabs;
                    Fi += Ft;
                }
            }

            /* ─ 流体力 (stub) ─ */
            ApplyHydrodynamicForce(i,Fi);

            /* ─ 端点へ反力返送 or 内部節点 → RH ─ */
            if(i==0 || i==Seg){
                unsigned off=(i==0)?(Seg-1)*3:(Seg-1)*3+3;
                RH.Put(GetFirstIndex()+off+0,-Fi(1));
                RH.Put(GetFirstIndex()+off+1,-Fi(2));
                RH.Put(GetFirstIndex()+off+2,-Fi(3));
            }
        }
        return RH;
    }

    /* --- ヤコビアン （線形海底反力を追加） ------------------------------ */
    VariableSubMatrixHandler&
    AssJac(VariableSubMatrixHandler& WH, doublereal dCoef,
           const VectorHandler&, const VectorHandler&) override
    {
        const unsigned dim=(Seg-1)*3;
        DenseSubMatrixHandler& K = WH.SetDense(dim,dim); K.Zero();

        for(unsigned i=1;i<Seg;i++){
            unsigned r=(i-1)*3;
            const Prop left = GetMaterial(i-1);
            const Prop right= GetMaterial(i);

            /* 質量項: dCoef*m*I */
            const doublereal m = GetMaterial(i).m;
            for(unsigned a=0;a<3;a++) K(r+a,r+a)+=dCoef*m;

            /* 伸び剛性 (軸方向のみ) */
            K(r  ,r  ) += (left.EA/left.l0 + right.EA/right.l0);
            K(r+1,r+1)+= (left.EA/left.l0 + right.EA/right.l0);
            K(r+2,r+2)+= (left.EA/left.l0 + right.EA/right.l0);
            if(i>1)     for(unsigned a=0;a<3;a++) K(r+a-3,r+a)-=left.EA/left.l0;
            if(i<Seg-1) for(unsigned a=0;a<3;a++) K(r+a+3,r+a)-=right.EA/right.l0;

            /* 海底ばね dFn/dz */
            const Vec3 pos = N[i]->GetXCurr();
            const doublereal zb = pBed->GetZ(pos(1),pos(2));
            if(pos(3)<zb) K(r+2,r+2)+=kn;
        }
        return WH;
    }

    /* --- 出力 ----------------------------------------------------------- */
    void Output(OutputHandler& OH) const override
    {
        for(unsigned i=0;i<=Seg;i++){
            OH << "LM " << i << ' ' << N[i]->GetXCurr() << '\n';
        }
    }
};

/* === Factory & 登録 ======================================================= */
static UserDefinedElem*
CreateModuleLM(unsigned uLabel,const DofOwner* pDO,
               DataManager* pDM, MBDynParser& HP)
{
    return new ModuleLM(uLabel,pDO,pDM,HP);
}

extern "C" int module_init()
{
    bool ok=true;
    ok &= SetElemCreators("lumped_mooring",CreateModuleLM);
    ok &= SetElemCreators("catenary",CreateModuleLM);      /* 旧名互換 */
    if(!ok){
        silent_cerr("module-lm: element name already used"<<std::endl);
        return -1;
    }
    return 0;
}
