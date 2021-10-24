/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**parse  mapping result SAM file to figure out how many reads mapped to known genes
 * and calculate RPKMs, plot
 * requires Bedtools to be installed
 *
 * @author joey
 */
public class StepMappedReadsCountingInGenes extends NGSStep{
    static Logger                       logger = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS               = NGSStepSubclass.DATABLE;
    
    public static final String          STEP_ID_STRING              = "GenesCountOfMappedReads";
    private static final String         ID_BED_TOOLS                = "bedtools";
    private static final String         ID_SAM_TOOLS                = "samtools";
    private static final String         ID_KNOWN_GENES_BED          = "knownGenesInBED";
    

    
    
    private static final String         INFILE_EXTENSION            = ".trim.gen.sam";
    private static final String         SORTED_SAM_EXTENSION            = ".trim.gen.sorted.sam";
    private static final String         BAM_EXTENSION            = ".trim.gen.bam";
    private static final String         OUTFILE_GENES_EXTENSION     = ".genes.bed.coverage";
    private static final String         OUTFILE_FORWARD_EXTENSION   = ".genes.bed.forward.coverage";
    private static final String         OUTFILE_REVERSE_EXTENSION   = ".genes.bed.reverse.coverage";
    
    private static final String        RAW_INPUT_EXTENSION          = ".fastq.gz";
    
    private             String          knownGenesBED              = "";
    private             String          bedtools              = "";
     private             String          samtools              = "";
   


    public StepMappedReadsCountingInGenes(){
        classSubtype = NGSStepSubclass.DATABLE;
    }
    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepMappedReadsCountingInGenes(InputDataForStep sid) {
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    

    @Override
    public String shortStepDescription(){
      return "summarize SAM file";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "summarize SAM file.\n"
              + "i.e., how many reads mapped to known genes,  calculate RPKMs, plot results\n"
              + "Requires BedTools to be installed."
              + "The location can be specified in the YAML configuration file\n";
    }



    
    @Override
    public void verifyInputData() throws IOException  {
        logger.info("verify input data");  

        this.setPaths();
        
       if (new File(this.getBEDtools()).exists() == false) {
            logger.error("bedtools is not found at location < " + this.getBEDtools() + ">");
            throw new IOException("bedtools is not found at location < " + this.getBEDtools() + ">");
        }
       
       if (new File(this.getSAMtools()).exists() == false) {
            logger.error("bedtools is not found at location < " + this.getSAMtools() + ">");
            throw new IOException("bedtools is not found at location < " + this.getSAMtools() + ">");
        }
       
       if (new File(this.getKnownGenes()).exists() == false) {
            logger.error("known genes in BED format file is not found at location < " + this.getKnownGenes() + ">");
            throw new IOException("mknown genes in BED format file is not found at location < " + this.getKnownGenes() + ">");
        }

      
    }

    @Override
    public void verifyOutputData() throws IOException {
        
    }

    @Override
    public void parseConfigurationData(HashMap configData) throws Exception {
        logger.info(STEP_ID_STRING + ": verify configuration data");

        if (configData.get(ID_BED_TOOLS) == null) {
            logger.error("<" + configData.get(ID_BED_TOOLS) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_BED_TOOLS) + "> : Missing Definition in Configuration File");
        }
        
        if (configData.get(ID_SAM_TOOLS) == null) {
            logger.error("<" + configData.get(ID_SAM_TOOLS) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_SAM_TOOLS) + "> : Missing Definition in Configuration File");
        }
        
        if (configData.get(ID_KNOWN_GENES_BED) == null) {
            logger.error("<" + ID_KNOWN_GENES_BED + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_KNOWN_GENES_BED + "> : Missing Definition in Configuration File");
        }
        
        this.setBEDTools((String) configData.get(ID_BED_TOOLS));
        this.setSAMTools((String) configData.get(ID_SAM_TOOLS));
        this.setKnownGenes((String) configData.get(ID_KNOWN_GENES_BED));

        logger.info("passed");
        

    }

    @Override
    public HashMap generateExampleConfigurationData() {
        logger.info(STEP_ID_STRING + ": generate example configuration data");
        
        HashMap configData = new HashMap();
        
        configData.put(ID_BED_TOOLS, "/usr/bin/bedtools");
        configData.put(ID_KNOWN_GENES_BED, "/home/joey/projects/forJavaDevlop/hgTables-19.bed");

        return configData;
    }

    @Override
    public void execute() throws IOException {
        logger.info(STEP_ID_STRING + ": execute");
        
        Boolean fA = new File(outFolder).mkdir();
        if (fA) {
            logger.info("created output folder <" + outFolder + "> for results");
        }
        
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()) {
            try {
                SampleDataEntry sampleData = (SampleDataEntry) itSD.next();
                sortSAMFile(sampleData);
                samTobam(sampleData);
                countReadsInGenes(sampleData);
            }
            catch(Exception exception){}
                logger.error("error executing counting mapped reads in known genes command\n");
                throw new IOException(STEP_ID_STRING + ": error executing counting mapped reads in known genes command ");
        } 
        logger.info(STEP_ID_STRING + ": completed");
    }

    @Override
    public NGSStepSubclass getStepSubclass() {
        return classSubtype;
    }

    private void setBEDTools(String bedtoolsString) {
        this.bedtools = bedtoolsString;
    }

    private void setKnownGenes(String knownGenesString) {
        this.knownGenesBED = knownGenesString;
    }

    private String getBEDtools() {
        return bedtools;
    }

    private String getKnownGenes() {
        return knownGenesBED;
    }

    private void sortSAMFile(SampleDataEntry sampleData) throws IOException {
        logger.info("sorting SAM file...");
        String sortSAMcmd = "";
        try {
            ArrayList<String> cmd = new ArrayList<>();
            cmd.add(this.getSAMtools());

            cmd.add("sort");
            cmd.add(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
            cmd.add("-o ");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, SORTED_SAM_EXTENSION));

            sortSAMcmd = this.cleanPath(StringUtils.join(cmd, " "));
            logger.info("SAM file sorting command:\t" + sortSAMcmd);

            Runtime rt = Runtime.getRuntime();
            Process proc = rt.exec(sortSAMcmd);
            BufferedReader brAStdin = new BufferedReader(new InputStreamReader(proc.getInputStream()));
            BufferedReader brAStdErr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));

            String line = null;
            logger.info("<OUTPUT>");
            while ((line = brAStdin.readLine()) != null) {
                logger.info(line);
            }
            logger.info("</OUTPUT>");
            
            logger.info("<ERROR>");
            while ((line = brAStdErr.readLine()) != null) {
                logger.info(line);   
            }
            logger.info("</ERROR>");

            int exitVal = proc.waitFor();
            logger.info("Process exitValue: " + exitVal);

            brAStdin.close();
            brAStdErr.close();
            logger.info("sorting SAM file done");
            
      
        } catch (IOException | InterruptedException ex) {
            logger.error("error while sorting SAM file\n");
            logger.error(sortSAMcmd);
            logger.error(ex.toString());
            throw new IOException(STEP_ID_STRING + ": \"error while sorting SAM file " + sortSAMcmd);
        }

    }

    private String getSAMtools() {
        return samtools;
    }

    private void setSAMTools(String samtoolsString) {
        this.samtools = samtoolsString;
    }

    private void countReadsInGenes(SampleDataEntry sampleData) throws IOException, InterruptedException {
        /**
         * since when use process.exec(), it could not recognize "|" in shell command
         * so split the command into several lines
         */
        logger.info("counting readss in known genes");
        String countCmdString = "";
        try {
            ArrayList<String> cmd = new ArrayList<>();
            /**
             * /usr/bin/bedtools bamtobed -i /home/joey/projects/forJavaDevlop/genesCountReads/testingReads01.trim.gen.bam 
             *  > /home/joey/projects/forJavaDevlop/genesCountReads/testingReads01.trim.gen.bed
             */
            cmd.add(this.getBEDtools());
            cmd.add("bamtobed -i");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, BAM_EXTENSION));
            cmd.add(">");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ".trim.gen.bed"));
            cmd.add(";");
            /**
             * /usr/bin/bedtools coverage -a /home/joey/projects/forJavaDevlop/genesCountReads/testingReads01.trim.gen.bed
             * -b /home/joey/projects/forJavaDevlop/hgTables-19.bed 
             * > /home/joey/projects/forJavaDevlop/genesCountReads/testingReads01.genes.bed.coverage
             */
            cmd.add(this.getBEDtools());
            cmd.add("coverage");
            cmd.add("-a");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ".trim.gen.bed"));
            cmd.add(" -b");
            cmd.add(knownGenesBED);
            cmd.add(">");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, OUTFILE_GENES_EXTENSION));
            
            cmd.add(";");
            /**
             * grep + /home/joey/projects/forJavaDevlop/genesCountReads/testingReads01.trim.gen.bed
             * > /home/joey/projects/forJavaDevlop/genesCountReads/testingReads01.trim.gen.forward.bed
             */
            cmd.add("grep + ");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ".trim.gen.bed"));
            cmd.add(">");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ".trim.gen.forward.bed"));
            
            cmd.add(";");
            /**
             * grep - /home/joey/projects/forJavaDevlop/genesCountReads/testingReads01.trim.gen.bed
             * > /home/joey/projects/forJavaDevlop/genesCountReads/testingReads01.trim.gen.reverse.bed
             */
            cmd.add("grep - ");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ".trim.gen.bed"));
            cmd.add(">");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ".trim.gen.reverse.bed"));
            
            cmd.add(";");
            /**
             * /usr/bin/bedtools coverage -a /home/joey/projects/forJavaDevlop/genesCountReads/testingReads01.trim.gen.forward.bed
             * -b /home/joey/projects/forJavaDevlop/hgTables-19.bed 
             * > /home/joey/projects/forJavaDevlop/genesCountReads/testingReads01.genes.bed.coverage
             */
            cmd.add(this.bedtools);
            cmd.add("coverage");
            cmd.add("-a");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ".trim.gen.forward.bed"));
            cmd.add(" -b");
            cmd.add(knownGenesBED);
            cmd.add(">");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, OUTFILE_FORWARD_EXTENSION));
            
            cmd.add(";");
            /**
             * /usr/bin/bedtools coverage -a /home/joey/projects/forJavaDevlop/genesCountReads/testingReads01.trim.gen.reverse.bed
             * -b /home/joey/projects/forJavaDevlop/hgTables-19.bed 
             * > /home/joey/projects/forJavaDevlop/genesCountReads/testingReads01.genes.bed.coverage
             */
            cmd.add(this.bedtools);
            cmd.add("coverage");
            cmd.add("-a");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ".trim.gen.reverse.bed"));
            cmd.add(" -b");
            cmd.add(knownGenesBED);
            cmd.add(">");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, OUTFILE_REVERSE_EXTENSION));
            

            countCmdString = this.cleanPath(StringUtils.join(cmd, " "));
            logger.info("counting readss command:\t" + countCmdString);

            Runtime rt = Runtime.getRuntime();
//            FileOutputStream outfile = new FileOutputStream("tmpShellCooamnd.sh");
//            outfile.write(countCmdString);
            FileUtils.writeStringToFile(new File((outFolder + FILESEPARATOR +"tmpShellCooamnd.sh")), countCmdString);
            Process proc = rt.exec(new String[] {"sh", (outFolder + FILESEPARATOR +"tmpShellCooamnd.sh")});
            BufferedReader brAStdin = new BufferedReader(new InputStreamReader(proc.getInputStream()));
            BufferedReader brAStdErr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));

            String line = null;
            logger.info("<OUTPUT>");
            while ((line = brAStdin.readLine()) != null) {
                logger.info(line);
            }
            logger.info("</OUTPUT>");
            
            logger.info("<ERROR>");
            while ((line = brAStdErr.readLine()) != null) {
                logger.info(line);   
            }
            logger.info("</ERROR>");

            int exitVal = proc.waitFor();
            logger.info("Process exitValue: " + exitVal);

            brAStdin.close();
            brAStdErr.close();
            logger.info("counting readss in known genes done");
        } catch (InterruptedException  ex) {
            logger.error("error while counting readss in known genes\n");
            logger.error(countCmdString);
            logger.error(ex.toString());
            throw new InterruptedException(STEP_ID_STRING + ": \"error while counting reads in known genes " + countCmdString);
        }
    }

    private void samTobam(SampleDataEntry sampleData) throws IOException {
              
            /**
             * samtools view -Sb  testingReads01.trim.gen.sorted.sam> test.bam
             */
            
            logger.info("converting sam file to bam file");
        String convertCmdString = "";
        try {
            ArrayList<String> cmd = new ArrayList<>();
            cmd.add(this.getSAMtools());
            cmd.add("view -Sb ");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, SORTED_SAM_EXTENSION));
            cmd.add("-o");
            cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, BAM_EXTENSION));

            convertCmdString = this.cleanPath(StringUtils.join(cmd, " "));
            logger.info("counting readss command:\t" + convertCmdString);

            Runtime rt = Runtime.getRuntime();
            Process proc = rt.exec(convertCmdString);
            BufferedReader brAStdin = new BufferedReader(new InputStreamReader(proc.getInputStream()));
            BufferedReader brAStdErr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));

            String line = null;
            logger.info("<OUTPUT>");
            while ((line = brAStdin.readLine()) != null) {
                logger.info(line);
            }
            logger.info("</OUTPUT>");
            
            logger.info("<ERROR>");
            while ((line = brAStdErr.readLine()) != null) {
                logger.info(line);   
            }
            logger.info("</ERROR>");

            int exitVal = proc.waitFor();
            logger.info("Process exitValue: " + exitVal);

            brAStdin.close();
            brAStdErr.close();
            logger.info("counting readss in known genes done");
        } catch (IOException | InterruptedException ex) {
            logger.error("error while counting reads in known genes\n");
            logger.error(convertCmdString);
            logger.error(ex.toString());
            throw new IOException(STEP_ID_STRING + ": \"error while counting reads in known genes " + convertCmdString);
        }
    }
    
}
