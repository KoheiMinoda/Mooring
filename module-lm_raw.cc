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

class ModuleLM : virtual public Elem, public UserDefinedElem {
    unsigned Seg;
    doublereal EA, CA;
    doublereal zBed, kn, cn, mu;
    std::vector<doublereal> L0;
    std::vector<doublereal> M;
    std::vector<StructNode *> N;
public:
    ModuleLM(unsigned uL, const DofOwner *pDO, DataManager *pDM, MBDynParser &HP) : Elem(uL, pDO), UserDefinedElem(uL) {
        unsigned fair = HP.GetInt();
        unsigned anchor = HP.GetInt();
        if (HP.IsKeyWord("seg")) Seg = HP.GetInt(); else Seg = 20;
        if (HP.IsKeyWord("EA")) EA = HP.GetReal(); else EA = 1.e8;
        if (HP.IsKeyWord("CA")) CA = HP.GetReal(); else CA = 5.e5;
        doublereal l0 = 45.;
        doublereal m0 = 120.;
        if (HP.IsKeyWord("l0")) l0 = HP.GetReal();
        if (HP.IsKeyWord("mass")) m0 = HP.GetReal();
        L0.assign(Seg, l0);
        M.assign(Seg - 1, m0);
        zBed = -2.e2;
        kn = 1.e7;
        cn = 1.e5;
        mu = 0.6;
        if (HP.IsKeyWord("seabed")) {
            zBed = HP.GetReal(); kn = HP.GetReal(); cn = HP.GetReal(); mu = HP.GetReal();
        }
        N.resize(Seg + 1);
        N[0] = dynamic_cast<StructNode *>(pDM->ReadNode(fair));
        N[Seg] = dynamic_cast<StructNode *>(pDM->ReadNode(anchor));
        for (unsigned i = 1; i < Seg; ++i) {
            Vec3 X0 = N[0]->GetXCurr() + (N[Seg]->GetXCurr() - N[0]->GetXCurr()) * (static_cast<doublereal>(i) / Seg);
            StructNode *pN = pDM->CreateStructNode(pDO, X0, Vec3(0., 0., 0.));
            N[i] = pN;
        }
    }
    virtual ~ModuleLM() {}
    virtual std::string GetElemType() const { return "ModuleMooringLM"; }
    virtual unsigned int iGetNumDof() const { return (Seg - 1) * 3; }
    virtual void WorkSpaceDim(integer *piNumRows, integer *piNumCols) const {
        *piNumRows = *piNumCols = (Seg - 1) * 3;
    }
    Vec3 Tension(unsigned i, unsigned j) const {
        Vec3 rij = N[j]->GetXCurr() - N[i]->GetXCurr();
        doublereal dist = rij.Norm();
        if (dist == 0.) return Vec3(0., 0., 0.);
        Vec3 t = rij / dist;
        doublereal dL = dist - L0[i];
        doublereal k = EA / L0[i];
        doublereal c = CA / L0[i];
        doublereal vrel = Dot(N[j]->GetVCurr() - N[i]->GetVCurr(), t);
        return t * (k * dL + c * vrel);
    }
    Vec3 SeabedForce(unsigned i) const {
        const Vec3 &x = N[i]->GetXCurr();
        const Vec3 &v = N[i]->GetVCurr();
        if (x(3) >= zBed) return Vec3(0., 0., 0.);
        doublereal dn = zBed - x(3);
        doublereal fn = kn * dn - cn * v(3);
        Vec3 fnVec(0., 0., fn);
        Vec3 vt = v - Vec3(0., 0., v(3));
        doublereal vtMag = vt.Norm();
        if (vtMag > 0.) {
            Vec3 ft = -vt / vtMag * std::min(mu * std::abs(fn), kt * vtMag);
            return fnVec + ft;
        }
        return fnVec;
    }
    virtual SubVectorHandler &AssRes(SubVectorHandler &RH, doublereal dCoef, const VectorHandler &, const VectorHandler &) {
        RH.SetSizeZero((Seg - 1) * 3);
        for (unsigned i = 1; i < Seg; ++i) {
            Vec3 F = Tension(i - 1, i) + Tension(i + 1, i) + Vec3(0., 0., -M[i - 1] * gVal) + SeabedForce(i) - M[i - 1] * N[i]->GetXPrimePrime();
            RH.Put(3 * (i - 1) + 1, F(1));
            RH.Put(3 * (i - 1) + 2, F(2));
            RH.Put(3 * (i - 1) + 3, F(3));
        }
        return RH;
    }
    virtual VariableSubMatrixHandler &AssJac(VariableSubMatrixHandler &WH, const VectorHandler &XCurr, const VectorHandler &XPrimeCurr) {
        const integer dim = (Seg - 1) * 3;
        DenseSubMatrixHandler &K = WH.SetDense(dim, dim);
        K.Zero();
        ColumnVector R0(dim);
        SubVectorHandler RH(dim);
        AssRes(RH, 0., XCurr, XPrimeCurr);
        for (integer i = 1; i <= dim; ++i) R0(i) = RH(i);
        for (integer k = 1; k <= dim; ++k) {
            VectorHandler Xpert = XCurr;
            Xpert.Put(k, Xpert(k) + epsFD);
            SubVectorHandler R1(dim);
            AssRes(R1, 0., Xpert, XPrimeCurr);
            for (integer r = 1; r <= dim; ++r) K(r, k) = (R1(r) - R0(r)) / epsFD;
            K(k, k) += 1.e3;
        }
        return WH;
    }
    virtual void Output(OutputHandler &OH) const {
        for (unsigned i = 0; i <= Seg; ++i) OH << GetLabel() << " " << i << " " << N[i]->GetXCurr() << std::endl;
    }
    virtual std::ostream &Restart(std::ostream &out) const { return out; }
    virtual doublereal dGetPrivData(unsigned int) const { return 0.; }
};

static UserDefinedElem *CreateModuleLM(unsigned uLabel, const DofOwner *pDO, DataManager *pDM, MBDynParser &HP) {
    return new ModuleLM(uLabel, pDO, pDM, HP);
}

static UserDefElemCreator regModuleLM("ModuleMooringLM", CreateModuleLM);
