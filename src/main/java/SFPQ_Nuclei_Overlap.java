import SFPQ_Nuclei_Overlap_Tools.Cell;
import SFPQ_Nuclei_Overlap_Tools.SFPQ_Nuclei_Overlap_Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.WaitForUserDialog;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.util.ImageProcessorReader;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import java.io.FileWriter;
import java.util.ArrayList;
import loci.plugins.in.ImporterOptions;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;



/**
 * Detect nuclei and SFPQ cells and compute their overlap
 * 
 * @author ORION-CIRB
 */
public class SFPQ_Nuclei_Overlap implements PlugIn {
    
    SFPQ_Nuclei_Overlap_Tools tools = new SFPQ_Nuclei_Overlap_Tools();
    private String imageDir = "";
    public  String outDirResults = "";
    public  String rootName = "";
    public BufferedWriter results;
   
    public void run(String arg) {
        try {
            if ((!tools.checkInstalledModules())) {
                return;
            } 

            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }
            
            // Find images with czi extension
            ArrayList<String> imageFile = tools.findImages(imageDir, "czi");
            if (imageFile == null) {
                IJ.showMessage("Error", "No images found with czi extension");
                return;
            }
            
            // Create output folder
            outDirResults = imageDir + File.separator + "Results" + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            // Write header for nuclei parameters file
            String header = "Image name\tNucleus ID\tNucleus area (µm2)\tSFPQ area (µm2)\tOverlap volume (µm2)\tOverlap/SFPQ volume (%)\n";
            FileWriter fwResults = new FileWriter(outDirResults + "results.xls", false);
            results = new BufferedWriter(fwResults);
            results.write(header);
            results.flush();
            
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFile.get(0));
            
            // Find image calibration
            tools.cal = tools.findImageCalib(meta);
            
            // Find channel names
            String[] channels = tools.findChannels(imageFile.get(0), meta, reader);

            // Channels dialog
            String[] chs = tools.dialog(channels);
            if (chs == null) {
                IJ.showStatus("Plugin canceled");
                return;
            }
            
            for (String f : imageFile) {
                rootName = FilenameUtils.getBaseName(f);
                System.out.println("--- ANALYZING IMAGE " + rootName + " ------");
                reader.setId(f);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);

                // Open nuclei channel
                System.out.println("- Analyzing " + chs[0] + " channel -");
                int indexCh = ArrayUtils.indexOf(channels, chs[0]);
                ImagePlus stackNucleus = BF.openImagePlus(options)[indexCh];
                ImagePlus imgNucleus = tools.doZProjection(stackNucleus, ZProjector.MAX_METHOD);
                tools.flush_close(stackNucleus); 
                
                // Find nuclei
                System.out.println("Finding nuclei....");
                Objects3DIntPopulation nucPop = tools.cellposeDetection(imgNucleus, true);
                System.out.println(nucPop.getNbObjects() + " nuclei found");
                
                // Open SFPQ channel
                System.out.println("- Analyzing " + chs[1] + " channel -");
                indexCh = ArrayUtils.indexOf(channels, chs[1]);
                ImagePlus stackSFPQ = BF.openImagePlus(options)[indexCh];
                ImagePlus imgSFPQ = tools.doZProjection(stackSFPQ, ZProjector.MAX_METHOD);
                tools.flush_close(stackSFPQ); 
                
                // Find SFPQ cells
                System.out.println("Finding SFPQ cells....");
                Objects3DIntPopulation sfpqPop = tools.cellposeDetection(imgSFPQ, false);
                System.out.println(sfpqPop.getNbObjects() + " SFPQ cells found");
                
                // Colocalize nuclei with SFPQ cells
                System.out.println("- Colocalizing nuclei with SFPQ cells -");
                ArrayList<Cell> cells = tools.findColocPop(nucPop, sfpqPop);
                
                // Draw results
                System.out.println("- Saving results -");
                tools.drawResults(cells, imgSFPQ, rootName, outDirResults);
                
                // Write results
                for(Cell cell: cells) {
                    results.write(rootName+"\t"+cell.params.get("label")+"\t"+cell.params.get("nucArea")+"\t"+cell.params.get("sfpqArea")+"\t"+
                            cell.params.get("overlapArea")+"\t"+cell.params.get("overlapPerc")+"\n");
                    results.flush();
                }
                tools.flush_close(imgNucleus); 
                tools.flush_close(imgSFPQ);
            }
            results.close();
        } catch (IOException | DependencyException | ServiceException | FormatException  ex) {
            Logger.getLogger(SFPQ_Nuclei_Overlap.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        System.out.println("--- All done! ---");
    }
}

           