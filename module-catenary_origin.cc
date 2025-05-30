/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "mbconfig.h"           // This goes first in every *.c,*.cc file

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <limits>

#include "dataman.h"
#include "userelem.h"
#include "module-catenary.h" // 自前のヘッダ

class ModuleCatenary
    // Elem クラスの仮想継承と UserDefinedElem クラスの通常継承．仮想継承はさらに別の祖先と共有されても一つにまとまる．通常継承は基底クラスの公開メンバーが派生クラスでも公開される
	: virtual public Elem, public UserDefinedElem {
		public:

			// ----コンストラクタとデストラクタ----
			// カテナリー要素を初期化するコンストラクタ
			ModuleCatenary(unsigned uLabel, const DofOwner *pDO, DataManager* pDM, MBDynParser& HP);

			// 仮想デストラクタ，リソースの解法を行う
			virtual ~ModuleCatenary(void);

			// ----シミュレーション関連メソッド----
			// シミュレーション結果を出力するためのメソッド
			virtual void Output(OutputHandler& OH) const;
			// ヤコビアン行列のサイズを計算
			virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

			// ヤコビアン行列を組み立てるメソッド
			VariableSubMatrixHandler& 
			AssJac(VariableSubMatrixHandler& WorkMat,
				doublereal dCoef, 
				const VectorHandler& XCurr,
				const VectorHandler& XPrimeCurr);

			// 残差ベクトルを組み立てるメソッド
			SubVectorHandler& 
			AssRes(SubVectorHandler& WorkVec,
				doublereal dCoef,
				const VectorHandler& XCurr, 
				const VectorHandler& XPrimeCurr);

			// ----データアクセスメソッド----
			// プライベートデータの数を返す
			unsigned int iGetNumPrivData(void) const;
			// 接続されているノードの数を返す
			int iGetNumConnectedNodes(void) const;
			// 接続されているノードのリストを取得
			void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
			// シミュレーション値を設定
			void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
				SimulationEntity::Hints *ph);
			// 再起動情報をストリームに出力
			std::ostream& Restart(std::ostream& out) const;

			// 初期化関連メソッド
			// 初期自由度の数を返す
			virtual unsigned int iGetInitialNumDof(void) const;
			// 初期ヤコビアン行列のサイズを計算
			virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
			// 初期ヤコビアン行列を組み立て
   			VariableSubMatrixHandler&
			InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		    	const VectorHandler& XCurr);
			// 初期残差ベクトルの組み立て
   			SubVectorHandler& 
			InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
			// ノード情報を取得するためのメソッド
    			void GetNode();

		private:
    		const StructNode	*g_pNode;	// フェアリーダーを表す MBDyn の構造ノードへのポインタ
			double GFPx,GFPy,GFPz;			// フェアリーダーの初期グローバル x,y,z 座標
			double APx, APy, APz;			// アンカーのグローバル x,y,z 座標
			double FPx, FPy, FPz;			// フェアリーダーのグローバル x,y,z 座標
			double L;						// ライン全長
			double S;						// ライン懸垂部長さ
			double w;						// チェーンの単重
			double L0_APFP;					// 水平張力が 0 となる時のアンカー・フェアリーダー間の水平距離
			double L_APFP;					// アンカー・フェアリーダー間の水平距離
			double H, V;					// 張力(水平方向・鉛直方向)
			double dFx, dFy, dFz;			// 張力(N)
			double Fx, Fy, Fz;				// 反力(N)
			double delta;					// 水平方向の変位量
			double x, d, l, h;				// パラメータ(funcd関数)
			double f, df, p0;				// パラメータ(funcd関数)
			double x1, x2, xacc;			// パラメータ(rtsafe関数)
			double X;						// パラメータ(逆双曲線関数)

			DriveOwner          FSF; // ランプアップに使う係数
			Vec3				F;			// 張力 F(x,y,z)
			Vec3				M;			// モーメント M(x,y,z)
			Vec3				GFP;		// フェアリーダーの初期座標
			Vec3				AP;			// アンカーの座標
			Vec3 				FP_node;	// フェアリーダーの座標
			Vec3				FP;			// フェアリーダーの座標
			Vec3				FP_AP;		// APを原点とした時のフェアリーダー座標
			Vec3				GFP_AP;		// APを原点とした時のフェアリーダー初期座標

		private:
			double myasinh(double X);
			double myacosh(double X);
			double myatanh(double X);

			void funcd(double x, double xacc, double &f, double &df, double d, double l, double &p0);
			double rtsafe(double x1, double x2, double xacc, double d, double l, double &p0);
	};

// クラス ModuleCatenary のコンストラクタをクラス定義の外側で実装している書き方
// ModuleCatenary:: はこの関数は ModuleCatenary クラスのメンバであるという合図
ModuleCatenary::ModuleCatenary(
	// コンストラクタに渡される引数リスト
	unsigned uLabel,
	const DofOwner *pDO, 
	DataManager* pDM, 
	MBDynParser& HP
)
	// コンストラクタ本体 { … } が実行される前にmここで列挙した「基底クラスのコンストラクタ呼び出し」や「メンバ変数の初期化」が行われる
	// Elem(uLabel, flag(0)) は ModuleCatenary が仮想継承している基底クラス Elem を呼び出すコンストラクタ, UserDefinedElem(uLabel, pDO) 通常継承しているもう一つの基底クラス UserDefinedElem のコンストラクタ
	: Elem(uLabel, flag(0)), UserDefinedElem(uLabel, pDO), 
		g_pNode(0), 
		GFPx(0), GFPy(0), GFPz(0), 
		APx(0), APy(0), APz(0),
		FPx(0), FPy(0), FPz(0),
		L(0), S(0), w(0),
		L0_APFP(0),	L_APFP(0),
		H(0), V(0),
		dFx(0), dFy(0), dFz(0),
		Fx(0), Fy(0), Fz(0),
		delta(0), x(0), d(0), l(0), h(0),
		p0(0), x1(0), x2(0), xacc(0), df(0), X(0)
	{
		// help
		if (HP.IsKeyWord("help")) {
			silent_cout(
				"\n"
				"Module: 	ModuleCatenary\n"
				"\n"
				<< std::endl
			);

			if (!HP.IsArg()) {
				/*
			 	* Exit quietly if nothing else is provided
			 	*/
				throw NoErr(MBDYN_EXCEPT_ARGS);
			}
		}

		// MBDyn のパーサー HP を介して各パラメータをメンバ変数に格納する
		g_pNode = dynamic_cast<const StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
		L  =  HP.GetReal();
		w  =  HP.GetReal();
		xacc =  HP.GetReal();
		APx =  HP.GetReal();
		APy =  HP.GetReal();
		APz =  HP.GetReal();

		// getting Force scale factor if needed.
		if (HP.IsKeyWord("Force" "scale" "factor")) {
			FSF.Set(HP.GetDriveCaller());

		} else {
			FSF.Set(new OneDriveCaller);
		}

		SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

		pDM->GetLogFile() << "catenary: "
			<< uLabel << " "
			<< std::endl;
	}

ModuleCatenary::~ModuleCatenary(void)
	{
		// destroy private data
		NO_OP;
	}

void ModuleCatenary::Output(OutputHandler& OH) const
	{
		if (bToBeOutput()) {

			if (OH.UseText(OutputHandler::LOADABLE)) {

				OH.Loadable() << GetLabel()
					<< " " << FPx
					<< " " << FPy
					<< " " << FPz
					<< " " << L_APFP
					<< " " << Fx
					<< " " << Fy
					<< " " << Fz
					<< std::endl;
			}
		}
	}

void ModuleCatenary::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const 
	{
		*piNumRows = 6;
		*piNumCols = 6;
	}

VariableSubMatrixHandler& ModuleCatenary::AssJac( VariableSubMatrixHandler& WorkMat, doublereal dCoef, const VectorHandler& XCurr, const VectorHandler& XPrimeCurr )
	{
		// should do something useful
		WorkMat.SetNullMatrix();

		return WorkMat;
	}

double ModuleCatenary::myasinh(double X) 
	{
		return std::log(X + std::sqrt(X * X + 1));
	}

double ModuleCatenary::myacosh(double X) 
	{
		return std::log(X + std::sqrt(X + 1) * std::sqrt(X - 1));
	}

double ModuleCatenary::myatanh(double X)
	{
		return 0.5 * std::log((1 + X) / (1 - X));
	}

// rtsafe 関数から呼び出され，カテナリー係留索の形状を決定するための方程式の残差とその導関数を計算する
// rtsafe 関数はこの funcd 関数が返す残差 f が 0 になるパラメータ x を見つけ出そうとする
// ここでの x (H/(wh)) は現時点での推定値
void ModuleCatenary::funcd(double x, double xacc, double& f, double& df, double d, double l, double& p0)
	{
    		int i,max;
		double f1, df1;
    		max = 1000;

		// 係留索が完全にたるんでいて，水平張力 0
    		if(x==0.0) {
        		f=-d;
        		df=0e-0;
        		p0=0e-0;
    		}

		// 水平張力あり
    		else if(x>0.0) {

			// 全長が垂直距離以下という物理的にあり得ない状況だが，特定の計算方法を適応する
        		if(l<=0.0) {
				double X_1;
				X_1 = 1.0/x+1.0;

				f=x*myacosh(X_1)-std::sqrt(1.0+2.0*x)+1.0-d;
				df=myacosh(X_1)-1.0/std::sqrt(1.0+2.0*x)-1.0/(x*std::sqrt(std::pow(X_1, 2.0)-1.0));
            			p0=0.0;
        		} else {
				// 海底に接する可能性のある複雑なケース
            			if(x>(l*l-1.0)/2) {
                			p0=0.0;
                			for(int i=1; i<max; i++) {
						double func1;
						func1 = 1.0/x+1.0/cos(p0);
					
						f1=x*(std::sqrt(std::pow(func1,2.0)-1.0)-std::tan(p0))-l;
						df1=x*(func1*std::tan(p0)*(1.0/cos(p0))/std::sqrt(std::pow(func1,2.0)-1.0)-std::pow(std::tan(p0), 2.0)-1.0);
                    				p0=p0-f1/df1;
						f1=x*(std::sqrt(std::pow(func1,2.0)-1.0)-std::tan(p0))-l;

                    				if(fabs(f1)<xacc) { break; }
                			}
				
					if(fabs(f1)>xacc) {
						std::cout<< "fabs(f1)>eps" << std::endl;
						throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
					}

					double X_2 = l/x+std::tan(p0);
					double X_3 = std::tan(p0);

					f=x*(myasinh(X_2)-myasinh(X_3))-l+1.0-d;
					df=myasinh(X_2)-myasinh(X_3)-l/(x*std::sqrt(std::pow(X_2, 2.0)+1.0));

				// 単純なカテナリー
				} else {
					double X_5;
					X_5 = 1.0/x+1.0;

					f=x*myacosh(X_5)-std::sqrt(1.0+2.0*x)+1.0-d;
					df=myacosh(X_5)-1.0/std::sqrt(1.0+2.0*x)-1.0/(x*std::sqrt(std::pow(X_5, 2.0)-1.0));
                			p0=0.0;
            			}
        		}
    		} else {
			std::cout << "ERROR (x<0)" << std::endl;
			throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
		}
	}

// カテナリー曲線の形状を決定するための方程式の解 H / (wh) を求めるために利用される
// x1 x2 : 初期区間の上限と下限 xacc : 要求される解の精度，これより小さくなれば収束判定 d l : カテナリー方程式を特徴づける無次元パラメータ &p0 : funcd 関数から返される角度パラメータ（参照渡し）
double ModuleCatenary::rtsafe(double x1, double x2, double xacc, double d, double l, double &p0)
	{
		const int MAXIT=1000;
    	int j;
		double fh,fl,xh,xl;
    	double dx,dxold,f,temp,rts;
		double p1, p2;

	// x1 x2 の各点で funcd を呼び出して解きたい方程式 f(x)=0 の f(x) とその導関数 df(x)/dx を計算
	// fl=f(x1), fh=f(x2) を代入
    	ModuleCatenary::funcd(x1, xacc, fl, df, d, l, p1);
    	ModuleCatenary::funcd(x2, xacc, fh, df, d, l, p2);

	// fl fh が同符号であれば根が無いため，中断する
    	if((fl>0.0&&fh>0.0)||(fl<0.0&&fh<0.0)) {
			std::cout << "ERROR (fl>0.0&&fh>0.0)||(fl<0.0&&fh<0.0)" << std::endl;
			throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
		}

	// fl または fh が 0 であれば x1 x2 が既に根であるため，その値を返す
    	if(fl==0.0) {
			p0  = p1;
			rts = x1;
			return rts;
		}
    	if(fh==0.0) {
			p0  = p2;
			rts = x2;
			return rts;
		}

	// f(xl)<0 かつ f(xh)>0 を満たすように根を含む区間の下限 xl と根を含む区間の上限 xh を設定
    	if(fl<0.0) {
        	xl=x1;
        	xh=x2;
    	} else {
        	xh=x1;
        	xl=x2;
    	}

	// 根の初期推定値は中点とする
    	rts=0.5*(x1+x2);
    	dxold=std::fabs(x2-x1);
    	dx=dxold;
	// funcd により初期推定値 rts における関数の値 f と微分 df を計算する
    	ModuleCatenary::funcd(rts, xacc, f, df, d, l, p0);

    	for(j=0; j<MAXIT; j++) {

		// ニュートン法のステップが区間の外に出るかどうか || ニュートン法のステップが二分法のステップよりも効果的かどうか
        	if((((rts-xh)*df-f)*((rts-xl)*df-f)>0.0)||((std::fabs(2.0*f))>std::fabs(dxold*df))) {
			// 二分法を採用した場合
            		dxold = dx;
            		dx = 0.5*(xh-xl);
            		rts =xl+dx;	
            		if(xl==rts) { return rts; }
        	} else {
			// ニュートン法を採用した場合
            		dxold=dx;
            		dx=f/df;
            		temp=rts;
            		rts-=dx;
            		if(temp==rts) {return rts;}
        	}

		// 収束判定
        	if(std::fabs(dx)<xacc) { return rts; }

		// 新しい推定値 rts で funcd を呼び出し，関数値 f を更新
		ModuleCatenary::funcd(rts, xacc, f, df, d, l, p0);

		// 区間の更新
		if(f<0.0){
            		xl=rts;
        	} else {
            		xh=rts;
        	}
    	}

		std::cout << "ERROR (Bisection method)" << std::endl;
		throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
	}

// 各タイムステップで呼び出され，カテナリー係留索要素がフェアリーダーノードに及ぼす力を計算する
// 運動方程式を解くことで次のタイムステップの状態を決定するが，AssRes 関数は外力に相当する部分を計算する
// WorkVec : この要素が計算する力とモーメントの寄与を格納するための MBDyn のシステム全体の残差ベクトルの一部を扱うハンドラ dCoef 時間積分スキームに関連する係数？ここでは使われていない？ XCurr システム全体の現在の変位を格納するベクトル XPrimeCurr システム全体の現在の速度を格納するベクトル
SubVectorHandler& ModuleCatenary::AssRes(SubVectorHandler& WorkVec, doublereal dCoef, const VectorHandler& XCurr, const VectorHandler& XPrimeCurr)
	{
		integer iNumRows = 0;					//ベクトルの大きさを決めるint変数を宣言して初期化
		integer iNumCols = 0;	
		WorkSpaceDim(&iNumRows, &iNumCols);		//先ほどのworkspecedim()関数で値を変更（6になった）
		WorkVec.ResizeReset(iNumRows);			//iNumRowsの大きさでWorkVecの大きさを変更し初期化

		/*
	 	*Fがどこに格納されるか示すインデックスをWorkVecに決める
	 	*/
		integer iFirstMomIndex = g_pNode->iGetFirstMomentumIndex();
		for (int iCnt = 1; iCnt<=6; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstMomIndex + iCnt);
		}

		const Vec3 FP = g_pNode->GetXCurr();			//ノード(浮体着鎖点)の現在の座標を取得		

		GFP = Vec3(GFPx,GFPy,GFPz);			// inputからGFPにフェアリーダーの初期座標を格納
		AP  = Vec3(APx,APy,APz);			// inputからAPにアンカーの初期座標を格納

		FPx = FP.dGet(1);
		FPy = FP.dGet(2);
		FPz = FP.dGet(3);

		// アンカーを原点とする座標変形
		FP_AP  = FP  - AP;
		GFP_AP = GFP - AP;

		// アンカー・フェアリーダー間の鉛直距離
		h = fabs(FP_AP.dGet(3));
	
		// 水平張力が0となる時のアンカー・フェアリーダー間の水平距離
		L0_APFP = L - h;
		// アンカー・フェアリーダー間の水平距離
		L_APFP  = std::sqrt(std::pow(FP_AP.dGet(1), 2)+std::pow(FP_AP.dGet(2), 2));

		// 関数に引き渡す値の計算
		delta = L_APFP-L0_APFP;
		d	  = delta / h;
		l	  = L / h;
	
		// 関数を用いた水平張力の計算
		if(d<=0) {
			//f=0;
			H  = 0;						// 水平張力
			V = w*h;					// 鉛直張力
			p0 = 0;
		} else if(d>=(std::sqrt(std::pow(l,2)-1)-(l-1))) {
			std::cout << "ERROR (The length between anchor to fairlead  is over the length of chain)" << std::endl;
			throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
		} else {		// d<(std::sqrt(std::pow(l,2)-1)-(l+1))
			x1 = 0;
			x2 = 1.e+6;
			//xacc = 1.e-6;
			// Ans_x = H/wh (H:水平張力, w:チェーンの単重, h:FPから海底までの距離)
			double Ans_x = ModuleCatenary::rtsafe(x1, x2, xacc, d, l, p0);

			S = h*std::sqrt(1+2*Ans_x);	// ライン懸垂部長さ
			H = Ans_x*w*h;				// 水平張力
			V = w*S;					// 鉛直張力
		}

		//drive callerによるramp up係数を乗じる
		doublereal dFSF = FSF.dGet();
		H = H*dFSF;
		V = V*dFSF;

		// ローカル座標系からグローバル座標系への変換
		if(FP_AP.dGet(1)>=0) {
			dFx = H*std::cos(std::atan((FP_AP.dGet(2))/(FP_AP.dGet(1))));
			dFy = H*std::sin(std::atan((FP_AP.dGet(2))/(FP_AP.dGet(1))));
			dFz = V;
		}
		// 固定着鎖点と浮体着鎖点のx座標が負の場合算出した角にπを足す
		else {
			dFx = H*std::cos(std::atan((FP_AP.dGet(2))/(FP_AP.dGet(1)))+M_PI);
			dFy = H*std::sin(std::atan((FP_AP.dGet(2))/(FP_AP.dGet(1)))+M_PI);
			dFz = V;
		}

		Fx = -dFx;
		Fy = -dFy;
		Fz = -dFz;

		F = Vec3(Fx,Fy,Fz);		// 張力 F(x,y,z)
		M = Vec3(0,0,0);			// モーメント　M(x,y,z)

		// 張力とモーメントをMBDynソルバーへ引き渡し
		WorkVec.Add(1,F);
		WorkVec.Add(4,M);

		return WorkVec;
	}

unsigned int ModuleCatenary::iGetNumPrivData(void) const
	{
		return 0;
	}

int ModuleCatenary::iGetNumConnectedNodes(void) const
	{
		return 0;
	}

void ModuleCatenary::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
	{
		connectedNodes.resize(0);
	}

void ModuleCatenary::SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP, SimulationEntity::Hints *ph)
	{
		NO_OP;
	}

std::ostream& ModuleCatenary::Restart(std::ostream& out) const
	{
		return out << "# ModuleTemplate: not implemented" << std::endl;
	}

unsigned int ModuleCatenary::iGetInitialNumDof(void) const
	{
		return 0;
	}

void ModuleCatenary::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
	{
		*piNumRows = 0;
		*piNumCols = 0;
	}

VariableSubMatrixHandler& ModuleCatenary::InitialAssJac( VariableSubMatrixHandler& WorkMat, const VectorHandler& XCurr)
	{
		// should not be called, since initial workspace is empty
		ASSERT(0);

		WorkMat.SetNullMatrix();

		return WorkMat;
	}

SubVectorHandler& ModuleCatenary::InitialAssRes( SubVectorHandler& WorkVec, const VectorHandler& XCurr)
	{
		// should not be called, since initial workspace is empty
		ASSERT(0);

		WorkVec.ResizeReset(0);

		return WorkVec;
	}

bool catenary_set(void) {
	#ifdef DEBUG
		std::cerr << __FILE__ <<":"<< __LINE__ << ":"<< __PRETTY_FUNCTION__ << std::endl;
	#endif

	UserDefinedElemRead *rf = new UDERead<ModuleCatenary>;

	if (!SetUDE("catenary", rf)) {
		delete rf;
		return false;
	}

	return true;
}

#ifndef STATIC_MODULES

extern "C" {

	int module_init(const char *module_name, void *pdm, void *php) {
		if (!catenary_set()) {
			silent_cerr("catenary: "
				"module_init(" << module_name << ") "
				"failed" << std::endl);
			return -1;
		}

		return 0;
	}

} // extern "C"

#endif // ! STATIC_MODULES
