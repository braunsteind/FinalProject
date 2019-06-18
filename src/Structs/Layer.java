package Structs;

public class Layer {
    public double[] dz;
    public double[] th_s;
    public double[] th_fc;
    public double[] th_wp;
    public double[] Ksat;
    public double[] Penetrability;
    public double[] th_dry;
    public double[] tau;

    public Layer() {
    }

    public Layer(Layer layer) {
        this.dz = layer.dz.clone();
        this.th_s = layer.th_s.clone();
        this.th_fc = layer.th_fc.clone();
        this.th_wp = layer.th_wp.clone();
        this.Ksat = layer.Ksat.clone();
        this.Penetrability = layer.Penetrability.clone();
        this.th_dry = layer.th_dry.clone();
        this.tau = layer.tau.clone();
    }
}