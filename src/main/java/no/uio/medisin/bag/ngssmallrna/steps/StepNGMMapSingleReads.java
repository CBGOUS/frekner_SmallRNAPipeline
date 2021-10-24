/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.ReferenceDataLocations;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import static no.uio.medisin.bag.ngssmallrna.steps.NGSStep.FILESEPARATOR;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * 
 * Use NextGenMap mapping tool for mapping reads that are likely to contain
 * high numbers of mismatches
 *
 * Input is a collapsed FASTQ file, not a FASTA file
 * 
 * @see <a href="http://cibiv.github.io/NextGenMap/">NextGenMap</a>
 * 
 * 
 * @author sr
 */


public class StepNGMMapSingleReads extends NGSStep{

    static Logger                       logger = LogManager.getLogger();

    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.DATABLE;

    public  static final String STEP_ID_STRING                  = "NGMMapSingleReads";

    private static final String ID_SOFTWARE                     = "mappingSoftware";
    private static final String ID_REF_SEQUENCE                 = "refSeq";
    private static final String ID_THREADS                      = "noOfThreads";
    private static final String ID_NO_OF_SCRIPTS                = "noOfScripts";

    private static final String INFILE_EXTENSION                = ".fastq";
    private static final String BAM_ALN_EXTENSION               = ".ngm.bam";
    private static final String BAM_SORT_EXTENSION              = ".ngm.sorted.bam";
    private static final String RAW_COUNTS_EXTENSION            = "raw.counts";
    private static final String MAPPED_COUNTS_EXTENSION         = ".mapped.counts";

    private             String  mappingSoftware                 = "";
    private             int     NoOfThreads                     = 2;
    private             int     NoOfScripts                     = 1;
    private             String  rootDataFolder                  = "";
    private             String  ReferenceGenome                 = "";
    private             String  pathToNGMRefIndex               = "";
    private             String  bamReferenceAln                 ="";
    
    private  ArrayList<String>  mapGenStdErr;
    private  ArrayList<String>  mapAbunStdErr;
    
    
    public StepNGMMapSingleReads() {
        classSubtype = NGSStepSubclass.DATABLE;
    }

    /**
     *
     * @param sid StepInputData
     *
     */
    public StepNGMMapSingleReads(InputDataForStep sid) {
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    
    
    @Override
    public String shortStepDescription(){
      return "Map reads with high numbers of mismatche to reference\n";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "Map reads with high numbers of mismatche to reference.\n"
              + "The step requires the NextGenMap software to be installed.\n"
              + "The location can be specified in the YAML configuration file \n"
              + "along with other parameters\n";
    }


    @Override
    public void parseStepParameters() throws Exception{
      
    }
    

    
    
    /**
     * in this method we are simply checking that the configuration file 
     * has all the entries we need. We dont check if the values are acceptable
     * that is the role of the NGSStep.
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{

        logger.info(STEP_ID_STRING + ": verify configuration data");
        
        if(configData.get(ID_SOFTWARE)==null) {
            logger.error("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_REF_SEQUENCE)==null) {
            logger.error("<" + ID_REF_SEQUENCE + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_REF_SEQUENCE + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_THREADS)==null) {
            logger.error("<" + ID_THREADS + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_THREADS + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_NO_OF_SCRIPTS)==null) {
            logger.error("<" + ID_NO_OF_SCRIPTS + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_NO_OF_SCRIPTS + "> : Missing Definition in Configuration File");
        }
        
        
        
        try{
            this.setNoOfThreads((Integer) configData.get(ID_THREADS));
        }
        catch(NumberFormatException exNm){
            logger.error(ID_THREADS + " <" + configData.get(ID_THREADS) + "> is not an integer");
            throw new NumberFormatException(ID_THREADS + " <" + configData.get(ID_THREADS) + "> is not an integer");
        }        
        if (this.getNoOfThreads() <= 0){
            logger.error(ID_THREADS + " <" + configData.get(ID_THREADS) + "> must be positive integer");
            throw new IllegalArgumentException(ID_THREADS + " <" + configData.get(ID_THREADS) + "> must be positive integer");
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


        this.setReferenceGenome((String) configData.get(ID_REF_SEQUENCE));
        if(this.getReferenceGenome().length() !=3 ){
            logger.error(ID_REF_SEQUENCE + " <" + configData.get(ID_REF_SEQUENCE) + "> must be a 3 letter string");            
            throw new IllegalArgumentException(ID_REF_SEQUENCE + " <" + configData.get(ID_REF_SEQUENCE) + "> must be a 3 letter string");            
        }
        this.setMappingSoftware((String) configData.get(ID_SOFTWARE));
        

        logger.info("passed");
    }
    
    
    
    
    
    /**
     * generate script to perform NGM mapping on specified files
     * basic form of commmand is:
     *   ngm 
     *   -r reference genome.fa (NGM will build its own index)
     *   -q query.fastq 
     *   -o align.sam  
     *   --bam 
     *   --threads
     * @throws IOException 
     */
    @Override
    public void execute() throws IOException{

        logger.info(STEP_ID_STRING + ": execute");        

        Boolean fA = new File(outFolder).mkdir();
        if (fA) {
            logger.info("created output folder <" + outFolder + "> for results");
        }
        logger.info("Mapping software is " + this.getMappingSoftware());
        
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
                            + ".NGMMAP." + cmdFileCount + ".sh");
                    BufferedWriter bwCmd = new BufferedWriter(new FileWriter(new File(cmdFileName)));
                        bwCmd.write(StringUtils.join(currentCmd, "\n"));

                    bwCmd.close();
                    
                    currentCmd.clear();
                    cmdFileCount++;
                    sampleCount=0;
                }
                
            } catch (IOException | InterruptedException ex) {
                logger.error("error writing mapping commands to <"+ cmdFileName + ">\n");
                logger.error(StringUtils.join(currentCmd, "\n"));
                logger.error(ex.toString());
                throw new IOException(STEP_ID_STRING + ": \"error writing mapping commands to < " + cmdFileName + ">");
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
        currentCmdString.add("#NGM map reads");
        currentCmdString.add(this.cmdMapReadsToReference(sampleData));
        currentCmdString.add("\n#sort BAM file");
        currentCmdString.add(this.cmdSortBAM(sampleData));
        currentCmdString.add("\n#count total reads");
        currentCmdString.add(this.cmdCountReadsInFastq(sampleData));
        currentCmdString.add("\n#count mapped reads");
        currentCmdString.add(this.cmdCountMappedReadsInBAM(sampleData));
        currentCmdString.add("\n\n");  
        
        return StringUtils.join(currentCmdString, "\n");
    }
    
    
    
    /**
     * summarize the mapping results for each sample
     * 
     * @return
     * @throws IOException 
     */
    private String summarizeMappingResults() throws IOException{
        
        
        DecimalFormat df = new DecimalFormat("###.##");

        String summaryString = "sample\traw counts\tmapped counts\tpercentage\n";
        String rawCountsFile = "";
        String mappedCountsFile = "";
        
        
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()) {
            SampleDataEntry sampleData = (SampleDataEntry) itSD.next();
            rawCountsFile = sampleData.getFastqFile1().replace(INFILE_EXTENSION, RAW_COUNTS_EXTENSION);  
            int rawCounts;
            int mappedCounts;
            
            try {
                
                BufferedReader brRC = new BufferedReader(new FileReader(new File(rawCountsFile)));
                    rawCounts = Integer.parseInt(brRC.readLine().trim());
                brRC.close();
                
                mappedCountsFile = sampleData.getFastqFile1().replace(INFILE_EXTENSION, MAPPED_COUNTS_EXTENSION);
                BufferedReader brMC = new BufferedReader(new FileReader(new File(mappedCountsFile)));
                    mappedCounts = Integer.parseInt(brRC.readLine().trim());
                brMC.close();
            } catch (IOException exIO) {
                logger.error("error opening result file <"+ rawCountsFile + ">\n");
                logger.error(StringUtils.join(rawCountsFile, "\n"));
                logger.error(exIO.toString());
                throw new IOException(STEP_ID_STRING + "error opening result file <"+ rawCountsFile + ">\n");
            }
            
    
            mappedCountsFile = sampleData.getFastqFile1().replace(INFILE_EXTENSION, MAPPED_COUNTS_EXTENSION);
            try {                
                BufferedReader brMC = new BufferedReader(new FileReader(new File(mappedCountsFile)));
                    mappedCounts = Integer.parseInt(brMC.readLine().trim());
                brMC.close();
            } catch (IOException exIO) {
                logger.error("error opening result file <"+ rawCountsFile + ">\n");
                logger.error(StringUtils.join(rawCountsFile, "\n"));
                logger.error(exIO.toString());
                throw new IOException(STEP_ID_STRING + "error opening result file <"+ rawCountsFile + ">\n");
            }
            
            summaryString = summaryString.concat(rawCounts 
                    + "\t" + mappedCounts 
                    + "\t" + df.format((double)mappedCounts/(double)rawCounts) + "\n");
        }
        
        return summaryString;
    }
        
        
        
        
    
    /**
     * map reads to the specified reference sequence
     * 
     * @param sampleData
     * @throws IOException
     * @throws InterruptedException 
     */
    private String cmdMapReadsToReference(SampleDataEntry sampleData) throws IOException, InterruptedException{
                 
        logger.info(STEP_ID_STRING + ": mapping reads to genome");
        String cmdNGMMapReads = "";
        ArrayList cmd = new ArrayList<>();
        cmd.add(this.getMappingSoftware());

        cmd.add("-r");
        cmd.add(pathToNGMRefIndex);

        String fastqInputFile = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1());
        cmd.add("-q " + fastqInputFile); 

        bamReferenceAln = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(INFILE_EXTENSION, BAM_ALN_EXTENSION));
        cmd.add("-o " + bamReferenceAln);

        cmd.add("--threads " + this.getNoOfThreads());

        cmd.add("-bam");



        cmdNGMMapReads = this.cleanPath(StringUtils.join(cmd, " "));
        logger.info("NGM Map Reads command:\t" + cmdNGMMapReads);

        return cmdNGMMapReads;


        /*
        
        try{
            Runtime rtGenMap = Runtime.getRuntime();
            Process procGenMap = rtGenMap.exec(cmdNGMMapReads);
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
            logger.error(cmdNGMMapReads);
            logger.error(ex.toString());
            throw new IOException(STEP_ID_STRING + ": \"error NGM Mapping genome reads " + cmdNGMMapReads);
        }

    */    
    }
    
    
    /**
     * count the lines in the fastq file associated with the sample entry
     * @param sampleData
     * @return line count
     * @throws IOException 
     */
    private long countReadsInFastq(SampleDataEntry sampleData) throws IOException{
        
        String fastqInputFile = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1());
        long lineCount = 0;
        try{
            LineNumberReader  lnr = new LineNumberReader(new FileReader(new File(fastqInputFile)));
            lnr.skip(Long.MAX_VALUE);
            lineCount= lnr.getLineNumber() + 1; 
            lnr.close(); 
        }
        catch(IOException exIO){
            logger.error("error counting lines in Fastq file\n");
            logger.error(exIO.toString());
            throw new IOException(STEP_ID_STRING + ": \"error counting lines in Fastq file <" + fastqInputFile + ">");
            
        }
        return lineCount/4;

    }
    
    
    /**
     * generate shell command to count number of lines in the fastq input file
     * 
     * @param sampleData
     * @return number of lines (not the number of reads)
     * @throws IOException 
     */
    private String cmdCountReadsInFastq(SampleDataEntry sampleData) throws IOException{
        return "wc -l " 
                + this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1() 
                + " > " 
                + sampleData.getFastqFile1().replace(INFILE_EXTENSION, RAW_COUNTS_EXTENSION)
    );
    }
    
    
    
    
    /**
     * generate shell commmand to sort the BAM file generated by the NGM mapping
     * 
     * @param sampleData
     * @return
     * @throws IOException
     * @throws InterruptedException 
     */
    private String cmdSortBAM(SampleDataEntry sampleData) throws IOException, InterruptedException{

        bamReferenceAln = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(INFILE_EXTENSION, BAM_ALN_EXTENSION));
        String sortedBAM = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(INFILE_EXTENSION, BAM_SORT_EXTENSION));
        String cmdBAMsort = "samtools sort " +  bamReferenceAln + " -o " + sortedBAM;
        return cmdBAMsort;
        /*
        Runtime rtGenMap = Runtime.getRuntime();
        Process procGenMap = rtGenMap.exec(cmdBAMsort);
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
     * count up the mapped reads in the generated BAM file
     * 
     * @param sampleData
     * @return numberOfMappedReads
     * @throws IOException
     * @throws InterruptedException 
     */
    private String cmdCountMappedReadsInBAM(SampleDataEntry sampleData) throws IOException, InterruptedException{

        bamReferenceAln = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(INFILE_EXTENSION, BAM_ALN_EXTENSION));
        String cmdBAMcount = "samtools view -F 0x904 -c " + bamReferenceAln + " > " 
                + sampleData.getFastqFile1().replace(INFILE_EXTENSION, MAPPED_COUNTS_EXTENSION);
        return cmdBAMcount;
        /*
        Runtime rtGenMap = Runtime.getRuntime();
        Process procGenMap = rtGenMap.exec(cmdBAMsort);
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
        pathToNGMRefIndex = this.cleanPath(getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_RIBOSOMAL_RNA);
        if(new File(pathToNGMRefIndex).exists() == false){
            logger.error("reference index not found at location < " + pathToNGMRefIndex +">");
            throw new IOException("reference index not found at location < " + pathToNGMRefIndex +">");
        }


        
        if(new File(this.getMappingSoftware()).exists() == false){
            logger.error("mapping software not found at location < " + this.getMappingSoftware() +">");
            throw new IOException("mapping software not found at location < " + this.getMappingSoftware() +">");
        }

        String pathToRefIndex = this.cleanPath(getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_RIBOSOMAL_RNA );
        if(new File(pathToRefIndex).isDirectory()== false){
            logger.error("abundant sequence bowtie index < " + pathToRefIndex +"> not found");
            throw new IOException("abundant sequence bowtie index  < " + pathToRefIndex +"> not found");
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
            
            // we dont check for Fastq 2 as this is single mapping
            
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

        configData.put(ID_SOFTWARE, "/home/user1/software/NextGenMap-0.5.0/bin/ngm-0.5.0/ngm");
        configData.put(ID_REF_SEQUENCE, "hsa");
        configData.put(ID_THREADS, 8);
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
     * @return the NoOfThreads
     */
    public int getNoOfThreads() {
        return NoOfThreads;
    }

    /**
     * @param NoOfThreads the NoOfThreads to set
     */
    public void setNoOfThreads(int NoOfThreads) {
        this.NoOfThreads = NoOfThreads;
    }

    /**
     * @return the ReferenceGenome
     */
    public String getReferenceGenome() {
        return ReferenceGenome;
    }

    /**
     * @param ReferenceGenome the ReferenceGenome to set
     */
    public void setReferenceGenome(String ReferenceGenome) {
        this.ReferenceGenome = ReferenceGenome;
    }

    /**
     * @return the mappingSoftware
     */
    public String getMappingSoftware() {
        return mappingSoftware;
    }

    /**
     * @param mappingSoftware the mappingSoftware to set
     */
    public void setMappingSoftware(String mappingSoftware) {
        this.mappingSoftware = mappingSoftware;
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


}

