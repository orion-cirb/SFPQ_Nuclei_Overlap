package SFPQ_Nuclei_Overlap_Tools;

import SFPQ_Nuclei_Overlap_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import SFPQ_Nuclei_Overlap_Tools.Cellpose.CellposeTaskSettings;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.Objects3DIntPopulationComputation;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.geom2.measurementsPopulation.PairObjects3DInt;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;


/**
 * @author ORION-CIRB
 */
public class SFPQ_Nuclei_Overlap_Tools {
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
    public Calibration cal = new Calibration();
    public float pixVol = 0;
    String[] chNames = {"Nuclei", "SFPQ cells"};

    // Cellpose
    public int cellPoseDiameter = 60;
    public String cellPoseModel = "cyto";
    public String cellPoseEnvDirPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose" : "/opt/miniconda3/envs/cellpose";
    public double cellposeStitchTh = 0.5;
    
    // Nuclei
    private double minNucVol = 200;
    private double maxNucVol = 2000;
    
    // SFPQ cells
    private double minSFPQInt = 20;
    
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.showMessage("Error", "3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
        
    /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in " + imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
  
    /**
     * Find image calibration
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
        return(cal);
    }
    
    
     /**
     * Find channels name
     * @param imageName
     * @return 
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        ArrayList<String> channels = new ArrayList<>();
        int chs = reader.getSizeC();
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelName(0, n).toString());
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelName(0, n).toString());
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelFluor(0, n).toString());
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelExcitationWavelength(0, n).value().toString());
                break; 
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelExcitationWavelength(0, n).value().toString());
                break;        
            default :
                for (int n = 0; n < chs; n++)
                    channels.add(Integer.toString(n));

        }
        return(channels.toArray(new String[channels.size()]));         
    }
    
    
    /**
     * Generate dialog box
     */
    public String[] dialog(String[] channels) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
          
        gd.addMessage("Channels", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.blue);
        int index = 0;
        for (String chName: chNames) {
            gd.addChoice(chName + ": ", channels, channels[index]);
            index++;
        }
        
        gd.addMessage("Nuclei detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min nucleus volume (µm3):", minNucVol);
        gd.addNumericField("Max nucleus volume (µm3):", maxNucVol); 
        
        gd.addMessage("SFPQ cells detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min intensity:", minSFPQInt); 
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY calibration (µm):", cal.pixelWidth);
        gd.addNumericField("Z calibration (µm):", cal.pixelDepth);
        gd.showDialog();
        
        String[] ch = new String[chNames.length];
        for (int i = 0; i < chNames.length; i++)
            ch[i] = gd.getNextChoice();
       
        minNucVol = (float) gd.getNextNumber();
        maxNucVol = (float) gd.getNextNumber();
        minSFPQInt = gd.getNextNumber();
        
        cal.pixelWidth = cal.pixelHeight = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();
        pixVol = (float) (cal.pixelWidth*cal.pixelHeight*cal.pixelDepth);
        
        if(gd.wasCanceled())
            ch = null;
                
        return(ch);
    }
    
    
    /**
     * Flush and close an image
     */
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Look for all 3D cells in a Z-stack: 
     * - apply CellPose in 2D slice by slice 
     * - let CellPose reconstruct cells in 3D using the stitch threshold parameters
     */
    public Objects3DIntPopulation cellposeDetection(ImagePlus imgIn, boolean nuclei) throws IOException{
        ImagePlus img = new Duplicator().run(imgIn);

        // Define CellPose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(cellPoseModel, 1, cellPoseDiameter, cellPoseEnvDirPath);
        settings.setStitchThreshold(cellposeStitchTh);
        settings.useGpu(true);
       
        // Run CellPose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, img);
        ImagePlus imgOut = cellpose.run();
        imgOut.setCalibration(cal);
       
        // Get cells as a population of objects
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgOut));
        System.out.println(pop.getNbObjects() + " CellPose detections");
       
        // Filter cells by size and intensity
        pop = zFilterPop(pop);
        if(nuclei) {
            pop = new Objects3DIntPopulationComputation​(pop).getFilterSize​(minNucVol/pixVol, maxNucVol/pixVol);
            System.out.println(pop.getNbObjects() + " detections remaining after size filtering");
        } else {
            filterDetectionsByIntensity(pop, img, minSFPQInt);
            System.out.println(pop.getNbObjects() + " detections remaining after intensity filtering");
        }
        pop.resetLabels();
        
        flush_close(img);
        flush_close(imgOut);
        return(pop);
    } 
    
    
    /*
     * Remove objects present in only one z slice from population 
     */
    public Objects3DIntPopulation zFilterPop (Objects3DIntPopulation pop) {
        Objects3DIntPopulation popZ = new Objects3DIntPopulation();
        for (Object3DInt obj : pop.getObjects3DInt()) {
            int zmin = obj.getBoundingBox().zmin;
            int zmax = obj.getBoundingBox().zmax;
            if (zmax != zmin)
                popZ.addObject(obj);
        }
        return popZ;
    }
    

    /**
     * Filter cells by intensity
     */
    public void filterDetectionsByIntensity(Objects3DIntPopulation cellPop, ImagePlus img, double intTh) {
        cellPop.getObjects3DInt().removeIf(p -> 
                (new MeasureIntensity(p, ImageHandler.wrap(img)).getValueMeasurement(MeasureIntensity.INTENSITY_AVG) < intTh));
    }
    
    
    /**
     * Find coloc objects in pop1 colocalized with pop2
     */
    public ArrayList<Cell> findColocPop (Objects3DIntPopulation pop1, Objects3DIntPopulation pop2, ImagePlus img) {
        AtomicInteger label = new AtomicInteger(0);
        ImageHandler imh = ImageHandler.wrap(img);
        ArrayList<Cell> cells = new ArrayList<>();
        if (pop1.getNbObjects() > 0 && pop2.getNbObjects() > 0) {
            MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(pop1, pop2);
            pop1.getObjects3DInt().forEach(obj1 -> {
                List<PairObjects3DInt> list = coloc.getPairsObject1(obj1.getLabel(), true);
                if (!list.isEmpty()) {
                    PairObjects3DInt p = list.get(list.size()-1);
                    double colocVol = p.getPairValue();
                    
                    if(colocVol > 0.5*obj1.size()) {
                        Object3DInt obj2 = p.getObject3D2();
                        Cell cell = new Cell(obj1, obj2);
                        
                        double sfpqInt = new MeasureIntensity(obj2, imh).getValueMeasurement(MeasureIntensity.INTENSITY_AVG);
                        cell.setParams(label.incrementAndGet(), obj1.size()*pixVol, obj2.size()*pixVol, sfpqInt, colocVol*pixVol);
                        
                        cells.add(cell);
                    }
                }
            });
        }
        return(cells);
    }
    
    
    /**
     * Save image objects
     */
    public void drawResults(ArrayList<Cell> cells, ImagePlus img, String imgName, String outDir) {
        ImageHandler imgObj1 = ImageHandler.wrap(new Duplicator().run(img)).createSameDimensions();
        ImageHandler imgObj2 = imgObj1.duplicate();
        
        for(Cell cell: cells) {
            float label = cell.params.get("label").floatValue();
            cell.nucleus.drawObject(imgObj1, label);
            cell.sfpq.drawObject(imgObj2, label);
        }
        
        ImagePlus[] imgColors = {imgObj2.getImagePlus(), null, imgObj1.getImagePlus(), img};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(cal);
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(outDir + imgName + ".tif"); 
        
        imgObj1.closeImagePlus();
        imgObj2.closeImagePlus();
        flush_close(imgObjects);
    }
    
}