package SFPQ_Nuclei_Overlap_Tools;

import java.util.HashMap;
import mcib3d.geom2.Object3DInt;


/**
 * @author ORION-CIRB
 */
public class Cell {
    
    public Object3DInt nucleus;
    public Object3DInt sfpq;
    public HashMap<String, Double> params;
    
    public Cell(Object3DInt nucleus, Object3DInt sfpq) {
        this.nucleus = nucleus;
        this.sfpq = sfpq;
        this.params = new HashMap<>();
    }
    
    public void setParams(double label, double nucVol, double sfpqVol, double sfpqInt, double overlapVol) {
        params.put("label", label);
        params.put("nucVol", nucVol);
        params.put("sfpqVol", sfpqVol);
        params.put("sfpqInt", sfpqInt);
        params.put("overlapVol", overlapVol);
        params.put("overlapPerc", overlapVol/sfpqVol*100);
    }
}