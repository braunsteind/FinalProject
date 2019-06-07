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
}
