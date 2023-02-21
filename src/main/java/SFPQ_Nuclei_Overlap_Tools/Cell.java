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
    
    public void setParams(double label, double nucArea, double sfpqArea, double overlapArea) {
        params.put("label", label);
        params.put("nucArea", nucArea);
        params.put("sfpqArea", sfpqArea);
        params.put("overlapArea", overlapArea);
        params.put("overlapPerc", overlapArea/sfpqArea*100);
    }
}