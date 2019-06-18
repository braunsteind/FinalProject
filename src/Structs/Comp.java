package Structs;

public class Comp {
    public Double[] dz;
    public Double[] dzsum;
    public Integer[] layer;
    public double[] th_fc;

    public Comp() {
    }

    public Comp(Comp comp) {
        this.dz = comp.dz.clone();
        this.dzsum = comp.dzsum.clone();
        this.layer = comp.layer.clone();
        this.th_fc = comp.th_fc.clone();
    }
}