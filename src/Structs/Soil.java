package Structs;

public class Soil {
    public double CalcSHP;
    public double Zsoil;
    public int nComp;
    public int nLayer;
    public double AdjREW;
    public double REW;
    public double CN;
    public double zRes;
    public double EvapZsurf;
    public double EvapZmin;
    public double EvapZmax;
    public double Kex;
    public double fevap;
    public double fWrelExp;
    public double fwcc;
    public double zCN;
    public double zGerm;
    public double AdjCN;
    public double fshape_cr;
    public double zTop;
    public Layer layer = new Layer();
    public Comp comp = new Comp();

    public Soil() {
    }

    ;

    public Soil(Soil soil) {
        this.CalcSHP = soil.CalcSHP;
        this.Zsoil = soil.Zsoil;
        this.nComp = soil.nComp;
        this.nLayer = soil.nLayer;
        this.AdjREW = soil.AdjREW;
        this.REW = soil.REW;
        this.CN = soil.CN;
        this.zRes = soil.zRes;
        this.EvapZsurf = soil.EvapZsurf;
        this.EvapZmin = soil.EvapZmin;
        this.EvapZmax = soil.EvapZmax;
        this.Kex = soil.Kex;
        this.fevap = soil.fevap;
        this.fWrelExp = soil.fWrelExp;
        this.fwcc = soil.fwcc;
        this.zCN = soil.zCN;
        this.zGerm = soil.zGerm;
        this.AdjCN = soil.AdjCN;
        this.fshape_cr = soil.fshape_cr;
        this.zTop = soil.zTop;
        this.layer = new Layer(soil.layer);
        this.comp = new Comp(soil.comp);
    }
}
