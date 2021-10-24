/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import static no.uio.medisin.bag.ngssmallrna.steps.NGSStep.FILESEPARATOR;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * 
 * Calculate site entropy from a BAM file using bamutils
 *
 * Input is a collapsed FASTQ file, not a FASTA file
 * 
 * @see <a href="http://cibiv.github.io/NextGenMap/">NextGenMap</a>
 * 
 * 
 * @author sr
 */


public class StepCalculateEntropyInBAM extends NGSStep{

    static Logger                       logger = LogManager.getLogger();

    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.DATABLE;

    public  static final String STEP_ID_STRING                  = "CalculateEntropyFromBAM";

    private static final String ID_BAMUTILS_SOFTWARE            = "bamutilsLocation";
    private static final String ID_PYTHON_SOFTWARE              = "pythonReformatter";
    private static final String ID_NO_OF_SCRIPTS                = "noOfScripts";
    private static final String ID_REF_SEQUENCE                 = "refSeq";

    private static final String INFILE_EXTENSION                = ".fastq";
    private static final String BAM_ALN_EXTENSION               = ".ngm.sorted.bam";
    private static final String BAM_INDEX_EXTENSION             = ".ngm.sorted.bam.bai";
    private static final String BAMUTILS_ENTROPY_EXTENSION      = ".bamutils.entropy.ngm.tsv";
    private static final String PROCESSED_ENTROPY_EXTENSION     = ".entropy.ngm.tsv";
    

    private             String  bamutilsSoftware                = "";
    private             String  pythonEntropyReformatter        = "";
    private             String  pathToChomosomeList             = "";
    private             String  rootDataFolder                  = "";
    private             int     NoOfScripts                     = 1;
    
    private             String  sortedBamReference                    ="";
    
    private  ArrayList<String>  mapGenStdErr;
    private  ArrayList<String>  mapAbunStdErr;
    
    
    public StepCalculateEntropyInBAM() {
        classSubtype = NGSStepSubclass.DATABLE;
    }

    /**
     *
     * @param sid StepInputData
     *
     */
    public StepCalculateEntropyInBAM(InputDataForStep sid) {
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    
    @Override
    public String shortStepDescription(){
      return "Calculate site entropy from a BAM file";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "Calculate site entropy from a BAM file.\n\n"
              + ""
              + "The step requires the bamutils software to be installed.\n"
              + "The location can be specified in the YAML configuration file\n"
              + "The step uses FASTQC in single mode, rather than batch analysis "
              + "so that a separate report is generated for each each file\n";
    }
  
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    
    /**
     * in this method we are simply checking that the configuration file 
     * has all the entries we need. We don't check if the values are acceptable,
     * that is the role of the NGSStep.
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{

        logger.info(STEP_ID_STRING + ": verify configuration data");
        

        if(configData.get(ID_BAMUTILS_SOFTWARE)==null) {
            logger.error("<" + ID_BAMUTILS_SOFTWARE + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_BAMUTILS_SOFTWARE + "> : Missing Definition in Configuration File");
        }
                
        this.setBamutilsSoftware((String) configData.get(ID_BAMUTILS_SOFTWARE));

        
        if(configData.get(ID_PYTHON_SOFTWARE)==null) {
            logger.error("<" + ID_PYTHON_SOFTWARE + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_PYTHON_SOFTWARE + "> : Missing Definition in Configuration File");
        }
                
        this.setPythonEntropyReformatter((String) configData.get(ID_PYTHON_SOFTWARE));
        
        if(configData.get(ID_NO_OF_SCRIPTS)==null) {
            logger.error("<" + ID_NO_OF_SCRIPTS + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_NO_OF_SCRIPTS + "> : Missing Definition in Configuration File");
        }

        
        try{
            this.setNoOfScripts((Integer) configData.get(ID_NO_OF_SCRIPTS));
        }
        catch(NumberFormatException exNm){
            logger.error(ID_NO_OF_SCRIPTS + " <" + configData.get(ID_NO_OF_SCRIPTS) + "> is not an integer");
            throw new NumberFormatException(ID_NO_OF_SCRIPTS + " <" + configData.get(ID_NO_OF_SCRIPTS) + "> is not an integer");
        }        
        if (this.getNoOfScripts() <= 0){
            logger.error(ID_NO_OF_SCRIPTS + " <" + configData.get(ID_NO_OF_SCRIPTS) + "> must be positive integer");
            throw new IllegalArgumentException(ID_NO_OF_SCRIPTS + " <" + configData.get(ID_NO_OF_SCRIPTS) + "> must be positive integer");
        }

        
        logger.info("passed");
    }
    
    
    
    
    
    /**
     * generate script to calculate entropy of mapped reads within a BAM file
     * 
     * basic form of commmand is:
     *      ./bamutils basecall 
     *      sample.sorted.bam 
     *      > sample.entropy.tsv
     * 
     * However, because of the structure of the code, it is necessary to switch
     * to the working directory for the software and execute from there.
     * 
     * @throws IOException 
     */
    @Override
    public void execute() throws IOException{

        logger.info(STEP_ID_STRING + ": execute");        

        Boolean fA = new File(outFolder).mkdir();
        if (fA) {
            logger.info("created output folder <" + outFolder + "> for results");
        }
        logger.info("Entropy software is " + this.getBamutilsSoftware());
        
        int sampleCount = 0;
        int cmdFileCount = 0;
        int mapsPerScript = this.getStepInputData().getSampleData().size()/this.getNoOfScripts();
        
        String cmdFileName = "";
        
        
        ArrayList<String> currentCmd = new ArrayList<>();
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()) {
            try {
                SampleDataEntry sampleData = (SampleDataEntry) itSD.next();
                sampleCount ++;
                currentCmd.add(buildCmdString(sampleData));
                
                if(sampleCount == mapsPerScript){
                    cmdFileName = this.cleanPath(outFolder 
                            + FILESEPARATOR 
                            + this.stepInputData.getProjectID()
                            + ".ENTCALC." + cmdFileCount + ".sh");
                    BufferedWriter bwCmd = new BufferedWriter(new FileWriter(new File(cmdFileName)));
                        bwCmd.write(StringUtils.join(currentCmd, "\n"));

                    bwCmd.close();
                    
                    currentCmd.clear();
                    cmdFileCount++;
                    sampleCount=0;
                }
                
            } catch (IOException | InterruptedException ex) {
                logger.error("error writing entropy commands to <"+ cmdFileName + ">\n");
                logger.error(StringUtils.join(currentCmd, "\n"));
                logger.error(ex.toString());
                throw new IOException(STEP_ID_STRING + ": \"error writing entropy commands to < " + cmdFileName + ">");
            }
            
        }
        if(sampleCount > 1){
            cmdFileName = this.cleanPath(outFolder 
                    + FILESEPARATOR 
                    + this.stepInputData.getProjectID()
                    + "." + cmdFileCount + ".sh");
            BufferedWriter bwCmd = new BufferedWriter(new FileWriter(new File(cmdFileName)));
                bwCmd.write(StringUtils.join(currentCmd, "\n"));

            bwCmd.close();

        }
        logger.info(STEP_ID_STRING + ": completed");
        
        

        /*
        it would be sensible here o summarize all the count information in a single file
        */
        
    }

    
    
    
    /**
     * create all the command strings for the current sample
     * 
     * @param sampleData
     * @return
     * @throws IOException
     * @throws InterruptedException 
     */
    private String buildCmdString(SampleDataEntry sampleData) throws IOException, InterruptedException{
 
        ArrayList<String> currentCmdString = new ArrayList<>();
        String sampleName = sampleData.getFastqFile1().replace(".fastq", "");
        currentCmdString.add("#+" + StringUtils.repeat("-", 60) + "+");                
        currentCmdString.add("#+" + StringUtils.repeat(" ", 60) + "+");
        currentCmdString.add("#+" + StringUtils.repeat(" ", 30-sampleName.length()/2) 
                + sampleName 
                + StringUtils.repeat(" ", 60 - (30-sampleName.length()/2) - sampleName.length())
                + "+" ) ;
        currentCmdString.add("#+" + StringUtils.repeat(" ", 60) + "+");
        currentCmdString.add("#+" + StringUtils.repeat("-", 60) + "+");

        currentCmdString.add("\n# samtools create index");
        currentCmdString.add("samtools index " 
                + this.cleanPath(inFolder + FILESEPARATOR 
                        + sampleData.getFastqFile1().replace(INFILE_EXTENSION, BAM_ALN_EXTENSION))
                + " " 
                + this.cleanPath(inFolder + FILESEPARATOR 
                        + sampleData.getFastqFile1().replace(INFILE_EXTENSION, BAM_INDEX_EXTENSION))
                );

        currentCmdString.add("# bamutils calc entropy");
        File f = new File(this.getBamutilsSoftware()); 
        currentCmdString.add("cd " + f.getParent());
        currentCmdString.add(this.cmdBamutilsCalcEntropy(sampleData));

        currentCmdString.add("\n# process bamutils entropy file");
        currentCmdString.add(this.cmdReformatEntropyResults(sampleData));

        currentCmdString.add("\n# generate summary plots in R");
        //currentCmdString.add(this.cmdPlotEntropyData(sampleData));
        currentCmdString.add("\n\n");  
        
        return StringUtils.join(currentCmdString, "\n");
    }
    
    /**
     * calculate entropy of mapped reads within the specified BAM file
     * 
     * basic form of commmand is:
     *      ./bamutils basecall 
     *      sample.sorted.bam 
     *      > sample.entropy.tsv
     *      @see <a href="http://ngsutils.org/modules/bamutils/basecall/"> http://ngsutils.org/modules/bamutils/basecall/</a>
     *      @param sampleData
     * @throws IOException
     * @throws InterruptedException 
     */
    private String cmdBamutilsCalcEntropy(SampleDataEntry sampleData) throws IOException, InterruptedException{
                 
        logger.info(STEP_ID_STRING + ": generate bamutils command to calculate entropy");
        String cmdBamutilsCalcEntropy = "";
        ArrayList cmd = new ArrayList<>();


        cmd.add("./bamutils basecall");

        cmd.add(this.cleanPath(inFolder 
                + FILESEPARATOR 
                + sampleData.getFastqFile1().replace(INFILE_EXTENSION, BAM_ALN_EXTENSION)));


        cmd.add(">");

        cmd.add(this.cleanPath(outFolder 
                + FILESEPARATOR 
                + sampleData.getFastqFile1().replace(INFILE_EXTENSION, BAMUTILS_ENTROPY_EXTENSION)));


        cmdBamutilsCalcEntropy = this.cleanPath(StringUtils.join(cmd, " "));
        logger.info("bamutils calc entropy command:\t" + cmdBamutilsCalcEntropy);

        return cmdBamutilsCalcEntropy;


        /*
        try{
            Runtime rtGenMap = Runtime.getRuntime();
            Process procGenMap = rtGenMap.exec(cmdBamutilsCalcEntropy);
            BufferedReader brGStdin = new BufferedReader(new InputStreamReader(procGenMap.getInputStream()));
            BufferedReader brGStdErr = new BufferedReader(new InputStreamReader(procGenMap.getErrorStream()));

            String line = "";
            String gLine = null;
            logger.info("<OUTPUT>");
            while ((gLine = brGStdin.readLine()) != null) {
                logger.info(line);
            }
            logger.info("</OUTPUT>");

            logger.info("<ERROR>");
            int errCount = 0;
            mapGenStdErr = new ArrayList<>();
            while ((line = brGStdErr.readLine()) != null) {
                logger.info(line);
                mapGenStdErr.add(line);
            }
            // need to parse the output from Bowtie to get the mapping summary
            logger.info(errCount + " errors reported");
            logger.info("</ERROR>");

            int gExitVal = procGenMap.waitFor();
            logger.info("Process exitValue: " + gExitVal);

            brGStdin.close();
            brGStdErr.close();
        } catch (IOException | InterruptedException ex) {
            logger.error("error NGM Mapping reads\n");
            logger.error(cmdBamutilsCalcEntropy);
            logger.error(ex.toString());
            throw new IOException(STEP_ID_STRING + ": \"error NGM Mapping genome reads " + cmdBamutilsCalcEntropy);
        }

    */    
    }
    
    
    /**
     * Reformat the data file generated by bamutils and recalculate the entropy
     * using a python script
     * 
     * Usage:
     * 
     *  process_entropy.py [-h] [-e ENTROPYFILELIST] [-c CHROMOSOMEFILE]
     * 
     * @param sampleData
     * @return number of lines (not the number of reads)
     * @throws IOException 
     */
    private String cmdReformatEntropyResults(SampleDataEntry sampleData) throws IOException{
        return this.getPythonEntropyReformatter()
                + " -c " + this.pathToChomosomeList
                + this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1());
    }
    
    
    
    /**
     * count up the mapped reads in the generated BAM file
     * 
     * @param sampleData
     * @return numberOfMappedReads
     * @throws IOException
     * @throws InterruptedException 
     */
    private String cmdPlotEntropyData(SampleDataEntry sampleData) throws IOException, InterruptedException{

        sortedBamReference = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(INFILE_EXTENSION, BAM_ALN_EXTENSION));
        String cmdBAMcount = "samtools view -F 0x904 -c " + sortedBamReference + " > " + sampleData.getFastqFile1().replace(INFILE_EXTENSION, ".counts.txt");
        return cmdBAMcount;
        /*
        Runtime rtGenMap = Runtime.getRuntime();
        Process procGenMap = rtGenMap.exec(cmdBAMcount);
        BufferedReader brGStdin = new BufferedReader(new InputStreamReader(procGenMap.getInputStream()));
        BufferedReader brGStdErr = new BufferedReader(new InputStreamReader(procGenMap.getErrorStream()));

        StringBuilder cmdOut = new StringBuilder();
        String line = "";
        String gLine = null;
        logger.info("<OUTPUT>");
        while ((gLine = brGStdin.readLine()) != null) {
            cmdOut.append(line);
        }
        logger.info("</OUTPUT>");

        logger.info("<ERROR>");
        int errCount = 0;
        mapGenStdErr = new ArrayList<>();
        while ((line = brGStdErr.readLine()) != null) {
            logger.info(line);
            mapGenStdErr.add(line);
        }
        
        
        return Integer.getInteger(cmdOut.toString());
        */
    }
    
    
    /**
     * @throws IOException
     * @throws NullPointerException 
     */
    @Override
    public void verifyInputData()  throws IOException, NullPointerException{

        logger.info(STEP_ID_STRING + ": verify input data");        
        this.setPaths();
        logger.info(STEP_ID_STRING + ": no path checking performed as "
            + "we are generating commands as shell scripts");
            
        /*
        pathToChomosomeList = this.cleanPath(getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() 
                + FILESEPARATOR + ReferenceDataLocations.ID_REL_RIBOSOMAL_RNA
                + FILESEPARATOR + ReferenceDataLocations.ID_GENOME_FEATURELIST_FILE);
        if(new File(pathToChomosomeList).exists() == false){
            logger.error("reference chromosome list not found at location < " + pathToChomosomeList +">");
            throw new IOException("reference chromosome list not found at location < " + pathToChomosomeList +">");
        }

        if(new File(this.getBamutilsSoftware()).exists() == false){
            logger.error("bamutils software not found at location < " + this.getBamutilsSoftware() +">");
            throw new IOException("bamutils software not found at location < " + this.getBamutilsSoftware() +">");
        }

        if(new File(this.getPythonEntropyReformatter()).exists() == false){
            logger.error("python software not found at location < " + this.getPythonEntropyReformatter() +">");
            throw new IOException("python software not found at location < " + this.getPythonEntropyReformatter() +">");
        }

        // does chromosome reference file exist?
        if(new File(this.getPythonEntropyReformatter()).exists() == false){
            logger.error("chromosome list reference was not found at  < " + this.getPythonEntropyReformatter() +">");
            throw new IOException("python software not found at location < " + this.getPythonEntropyReformatter() +">");
        }
        
                    
        // check the data files
        String fastqFile1;
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            
            //Fastq 1
            if (sampleData.getFastqFile1()==null) {
                logger.error(STEP_ID_STRING + ":no Fastq1 file specified");
                throw new IOException(STEP_ID_STRING + ":no Fastq1 file specified");
            }
            fastqFile1 = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(BAM_ALN_EXTENSION, INFILE_EXTENSION));
                        
            
            if (new File(this.cleanPath(fastqFile1)).exists()==false){
                logger.error(STEP_ID_STRING + ": fastq File1 <" 
                  + fastqFile1 + "> does not exist");
                throw new IOException(STEP_ID_STRING + ": fastq File1 <" 
                  + fastqFile1 + "> does not exist");
            }
            if (fastqFile1.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.error(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
            
            
        }
        logger.info("passed");        
        */
    }

    
    
    
    
    
    /**
     * generate sample configuration data so the user can see what can be
     * specified
     *
     * @return
     */
    @Override
    public HashMap generateExampleConfigurationData() {

        logger.info(STEP_ID_STRING + ": generate example configuration data");

        HashMap configData = new HashMap();

        configData.put(ID_BAMUTILS_SOFTWARE, "/Users/user1/software/ngsutils/bin/bamutils");
        configData.put(ID_PYTHON_SOFTWARE, "/Users/user1/software/python/bamutils_entropy_reformat.py");
        configData.put(ID_NO_OF_SCRIPTS, 2);

        return configData;
    }

    
    
    
    
    
    @Override
    public void verifyOutputData() {

    }

    @Override
    public NGSStepSubclass getStepSubclass(){
        return STEP_SUBCLASS;
    }
    
    
    /**
     * @return the bamutilsSoftware
     */
    public String getBamutilsSoftware() {
        return bamutilsSoftware;
    }

    /**
     * @return the rootDataFolder
     */
    public String getRootDataFolder() {
        return rootDataFolder;
    }

    /**
     * @param rootDataFolder the rootDataFolder to set
     */
    public void setRootDataFolder(String rootDataFolder) {
        this.rootDataFolder = rootDataFolder;
    }

    /**
     * @return the NoOfScripts
     */
    public int getNoOfScripts() {
        return NoOfScripts;
    }

    /**
     * @param NoOfScripts the NoOfScripts to set
     */
    public void setNoOfScripts(int NoOfScripts) {
        this.NoOfScripts = NoOfScripts;
    }

    /**
     * @param software the bamutilsSoftware to set
     */
    public void setBamutilsSoftware(String software) {
        this.bamutilsSoftware = software;
    }

    /**
     * @return the pythonEntropyReformatter
     */
    public String getPythonEntropyReformatter() {
        return pythonEntropyReformatter;
    }

    /**
     * @param pythonEntropyReformatter the pythonEntropyReformatter to set
     */
    public void setPythonEntropyReformatter(String pythonEntropyReformatter) {
        this.pythonEntropyReformatter = pythonEntropyReformatter;
    }

    /**
     * @return the pathToChomosomeList
     */
    public String getPathToChomosomeList() {
        return pathToChomosomeList;
    }

    /**
     * @param pathToChomosomeList the pathToChomosomeList to set
     */
    public void setPathToChomosomeList(String pathToChomosomeList) {
        this.pathToChomosomeList = pathToChomosomeList;
    }



}

