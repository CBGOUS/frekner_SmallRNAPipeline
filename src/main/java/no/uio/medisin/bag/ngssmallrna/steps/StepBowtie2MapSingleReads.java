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
import java.io.InputStreamReader;
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
 * Map reads to contaminants, RNA, Reference Genome and output in SAM format using bowtie2
 *
 * Input is a collapsed FASTA file
 * 
 * modified from class for bowtie mapping written by sr
 * 
 * @author zxf
 */


public class StepBowtie2MapSingleReads extends NGSStep {

    static Logger                       logger = LogManager.getLogger();

    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.DATABLE;

    public  static final String STEP_ID_STRING                  = "Bowtie2MapSingleReads";

    private static final String ID_SOFTWARE                     = "mappingSoftware";
    private static final String ID_REF_GENOME                   = "host";
    private static final String ID_MISMATCHES                   = "noOfMismatches";
    private static final String ID_ALIGN_MODE                   = "alignmentMode";
    private static final String ID_THREADS                      = "noOfThreads";

    private static final String INFILE_EXTENSION                = ".trim.clp.fasta";
    private static final String FASTQ_ABUNALN_EXTENSION         = ".trim.clp.abun.fasta";
    private static final String FASTQ_ABUNUNALN_EXTENSION       = ".trim.clp.notabun.fasta";
    private static final String SAM_ABUNALN_EXTENSION           = ".trim.clp.abun.sam";
    private static final String FASTQ_GENALN_EXTENSION          = ".trim.clp.gen.fasta";
    private static final String FASTQ_UNALN_EXTENSION           = ".trim.clp.unmap.fasta";
    private static final String SAM_GENALN_EXTENSION            = ".trim.clp.gen.sam";
    private static final String MAPPING_SUMMARY_EXTENSION       = ".trim.clp.gen.mapping.txt";
    
    private static final String        RAW_INPUT_EXTENSION      = ".fastq.gz";

    private             String  mappingSoftware                 = "";
    private             String  AlignMode                       = "";
    private             int     NoOfMismatches                  = 2;
    private             int     NoOfThreads                     = 4;
    private             String  rootDataFolder                  = "";
    private             String  ReferenceGenome                 = "";
    
    private             String  fastqTrimmedInputFile           = "";
    private             String  fastqAbundantAln                = "";
    private             String  fastqAbundantUnAln              = "";
    private             String  fastqGenomeAln                  = "";
    private             String  fastqGenomeUnAln                = "";
    private             String  samGenomeAln                    ="";
    
    private  ArrayList<String>  mapGenStdErr;
    private  ArrayList<String>  mapAbunStdErr;
    
    
    public StepBowtie2MapSingleReads() {
        classSubtype = NGSStepSubclass.DATABLE;
    }

    /**
     *
     * @param sid StepInputData
     *
     */
    public StepBowtie2MapSingleReads(InputDataForStep sid) {
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    

    @Override
    public String shortStepDescription(){
      return "Performs Bowtie2 mapping of single reads to reference genome";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "Performs Bowtie2 mapping of single reads to reference genome.\n"
              + "The step requires the Bowtie2 software to be installed.\n"
              + "The location can be specified in the YAML configuration file\n"
              + "along with mapping parameters and the reference genome";
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
        if(configData.get(ID_REF_GENOME)==null) {
            logger.error("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MISMATCHES)==null) {
            logger.error("<" + ID_MISMATCHES + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_MISMATCHES + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_ALIGN_MODE)==null) {
            logger.error("<" + configData.get(ID_ALIGN_MODE) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_ALIGN_MODE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_THREADS)==null) {
            logger.error("<" + ID_THREADS + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_THREADS + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_REF_GENOME)==null) {
            logger.error("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
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
            this.setNoOfMismatches((Integer) configData.get(ID_MISMATCHES));
        }
        catch(NumberFormatException exNm){
            logger.error(ID_MISMATCHES + " <" + ID_MISMATCHES + "> is not an integer");
            throw new NumberFormatException(ID_MISMATCHES + " <" + ID_MISMATCHES + "> is not an integer");
        }        
        if (this.getNoOfMismatches() <= 0){
            logger.error(ID_MISMATCHES + " <" + configData.get(ID_MISMATCHES) + "> must be positive integer");
            throw new IllegalArgumentException(ID_MISMATCHES + " <" + configData.get(ID_MISMATCHES) + "> must be positive integer");
        }
        
        

        this.setReferenceGenome((String) configData.get(ID_REF_GENOME));
        if(this.getReferenceGenome().length() !=3 ){
            logger.error(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
            throw new IllegalArgumentException(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
        }
        this.setAlignMode((String) configData.get(ID_ALIGN_MODE));
        this.setMappingSoftware((String) configData.get(ID_SOFTWARE));
        

        logger.info("passed");
    }
    
    
    
    
    
    /**
     * perform bowtie2 mapping on specified files
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
        String mappingCmd = this.getMappingSoftware();
        logger.info("Mapping software is " + mappingCmd);
        
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()) {
            ArrayList<String> cmd = new ArrayList<>();
            try {
                SampleDataEntry sampleData = (SampleDataEntry) itSD.next();

                this.mapAbundantReads(sampleData);
                this.mapReadsToGenome(sampleData);

                /*
                    write out mapping summary
                */
                String mappingOutputFile = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, MAPPING_SUMMARY_EXTENSION);
                BufferedWriter bwMO = new BufferedWriter(new FileWriter(new File(mappingOutputFile)));
                for (String mapLine : mapGenStdErr) {
                    bwMO.write(mapLine + "\n");
                }
                bwMO.write("\n\n" + "+" + StringUtils.repeat("-", 60) + "+" + "\n");
                bwMO.write("original FASTQ source" + sampleData.getFastqFile1() + "\n");
                bwMO.write(fastqTrimmedInputFile + "\n");
                bwMO.write(fastqGenomeAln + "\n");
                bwMO.write(fastqAbundantAln + "\n");
                bwMO.write(fastqGenomeUnAln + "\n");

                // Input 
                int totalInputReads = 0;
                String faLine = "";
                BufferedReader brIR = new BufferedReader(new FileReader(new File(fastqTrimmedInputFile)));
                while ((faLine = brIR.readLine()) != null) {
                    totalInputReads += Integer.parseInt(faLine.substring(1).split("-")[1]);
                    brIR.readLine();
                }
                brIR.close();
                bwMO.write("total input reads = " + totalInputReads + "\n");

                // Mapped
                int totalMappedReads = 0;
                faLine = "";
                BufferedReader brMR = new BufferedReader(new FileReader(new File(fastqGenomeAln)));
                while ((faLine = brMR.readLine()) != null) {
                    totalMappedReads += Integer.parseInt(faLine.substring(1).split("-")[1]);
                    brMR.readLine();
                }
                brMR.close();
                bwMO.write("total mapped reads = " + totalMappedReads + "\n");

                // Abundant
                int totalAbundantReads = 0;
                BufferedReader brAR = new BufferedReader(new FileReader(new File(fastqAbundantAln)));
                while ((faLine = brAR.readLine()) != null) {
                    totalAbundantReads += Integer.parseInt(faLine.substring(1).split("-")[1]);
                    brAR.readLine();
                }
                brAR.close();
                bwMO.write("total abundant reads = " + totalAbundantReads + "\n");

                // Unmapped
                int totalUnmappedReads = 0;
                BufferedReader brUR = new BufferedReader(new FileReader(new File(fastqGenomeUnAln)));
                while ((faLine = brUR.readLine()) != null) {
                    totalUnmappedReads += Integer.parseInt(faLine.substring(1).split("-")[1]);
                    brUR.readLine();
                }
                brUR.close();
                bwMO.write("total unmapped reads = " + totalUnmappedReads + "\n");

//                    bwMO.write("length filtered reads = " 
//                            + (rawReadsIn - totalMappedReads - totalAbundantReads - totalUnmappedReads));
                bwMO.write("\n\n" + "+" + StringUtils.repeat("-", 60) + "+" + "\n");

                bwMO.close();
            } catch (IOException | InterruptedException ex) {
                logger.error("error executing Bowtie Mapping command\n");
                logger.error(cmd);
                logger.error(ex.toString());
                throw new IOException(STEP_ID_STRING + ": \"error executing Bowtie Mapping command " + cmd);
            }
        }
        logger.info(STEP_ID_STRING + ": completed");

    }

    
    
    /**
     * Maps input reads to the supplied reference abundant sequences
     * 
     * @param sampleData
     * @throws IOException
     * @throws InterruptedException 
     */
    private void mapAbundantReads(SampleDataEntry sampleData) throws IOException, InterruptedException{
        
        logger.info(STEP_ID_STRING + ": mapping abundant reads");
        String cmdBowtieMapAbunReads = "";        
        try{
            ArrayList<String> cmd = new ArrayList<>();

            String mappingCmd = this.getMappingSoftware();


            cmd.add(mappingCmd);
            String pathToBowtieIndex = this.cleanPath(getStepInputData().getDataLocations().getGenomeRootFolder()
                    + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_ABUN_BOWTIE2_DATA_PATH + "abundantIndex2-all");
            cmd.add("-x");
            cmd.add(pathToBowtieIndex);

            fastqTrimmedInputFile = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
            cmd.add("-f");
            cmd.add(fastqTrimmedInputFile);

//            cmd.add("-" + this.getAlignMode());
//            cmd.add("--best");
//            cmd.add("-m " + this.getNoOfMismatches());

            fastqAbundantAln = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, FASTQ_ABUNALN_EXTENSION));
            fastqAbundantUnAln = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, FASTQ_ABUNUNALN_EXTENSION));
            String samAbundantAln = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, SAM_ABUNALN_EXTENSION));
            cmd.add("--al " + fastqAbundantAln);
            cmd.add("--un " + fastqAbundantUnAln);
            cmd.add("-p " + this.getNoOfThreads());
            cmd.add("");
            cmd.add("-S " + samAbundantAln);
            

            cmdBowtieMapAbunReads = this.cleanPath(StringUtils.join(cmd, " "));
            logger.info("Bowtie Map Abundant Reads command:\t" + cmdBowtieMapAbunReads);

            Runtime rt = Runtime.getRuntime();
            Process proc = rt.exec(cmdBowtieMapAbunReads);
            BufferedReader brAStdin = new BufferedReader(new InputStreamReader(proc.getInputStream()));
            BufferedReader brAStdErr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));

            String line = null;
            logger.info("<OUTPUT>");
            while ((line = brAStdin.readLine()) != null) {
                logger.info(line);
            }
            logger.info("</OUTPUT>");

            logger.info("<ERROR>");
            int skipCount = 0;
            mapAbunStdErr = new ArrayList<>();
            while ((line = brAStdErr.readLine()) != null) {
                if (line.contains("Warning: Skipping") && line.contains("less than")) {
                    skipCount++;
                } else {
                    logger.info(line);
                    mapAbunStdErr.add(line);
                }
            }

            // need to parse the output from Bowtie to get the mapping summary
            logger.info(skipCount + " lines were skipped because the read was too short");
            logger.info("</ERROR>");

            int exitVal = proc.waitFor();
            logger.info("Process exitValue: " + exitVal);

            brAStdin.close();
            brAStdErr.close();
            logger.info(STEP_ID_STRING + ": done");
        } catch (IOException | InterruptedException ex) {
            logger.error("error Bowtie Mapping abundant reads\n");
            logger.error(cmdBowtieMapAbunReads);
            logger.error(ex.toString());
            throw new IOException(STEP_ID_STRING + ": \"error Bowtie Mapping abundant reads " + cmdBowtieMapAbunReads);
        }
        
    }
    
    
    /**
     * map reads that didnt map to Abundant query sequences to the specified reference genome
     * 
     * @param sampleData
     * @throws IOException
     * @throws InterruptedException 
     */
    private void mapReadsToGenome(SampleDataEntry sampleData) throws IOException, InterruptedException{
                 
        logger.info(STEP_ID_STRING + ": mapping reads to genome");
        String cmdBowtieMapGenomeReads = "";
        try{
            String pathToBowtieGenomeIndex = this.cleanPath(getStepInputData().getDataLocations().getGenomeRootFolder()
                    + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_BOWTIE2_PATH);

            ArrayList cmd = new ArrayList<>();
            cmd.add(this.getMappingSoftware());
            cmd.add("-x");
            cmd.add(pathToBowtieGenomeIndex);

            cmd.add("-f");
            cmd.add(fastqAbundantUnAln);

//            cmd.add("-" + this.getAlignMode());
//            cmd.add("--best");
            cmd.add("-N " + this.getNoOfMismatches());
            
            cmd.add("-p " + this.getNoOfThreads());

            fastqGenomeAln = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, FASTQ_GENALN_EXTENSION));
            fastqGenomeUnAln = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, FASTQ_UNALN_EXTENSION));
//            fastqGenomeAln = this.cleanPath(outFolder + FILESEPARATOR);
//            fastqGenomeUnAln = this.cleanPath(outFolder + FILESEPARATOR );
            samGenomeAln = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, SAM_GENALN_EXTENSION));
            cmd.add("--al " + fastqGenomeAln);
            cmd.add("--un " + fastqGenomeUnAln);
            
            cmd.add("-S " + samGenomeAln);
            

            cmdBowtieMapGenomeReads = this.cleanPath(StringUtils.join(cmd, " "));
            logger.info("Bowtie Map Genome Reads command:\t" + cmdBowtieMapGenomeReads);

            Runtime rtGenMap = Runtime.getRuntime();
            Process procGenMap = rtGenMap.exec(cmdBowtieMapGenomeReads);
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
            int skipCount = 0;
            mapGenStdErr = new ArrayList<>();
            while ((line = brGStdErr.readLine()) != null) {
                if (line.contains("Warning: Skipping") && line.contains("less than")) {
                    skipCount++;
                } else {
                    logger.info(line);
                    mapGenStdErr.add(line);
                }
            }
            // need to parse the output from Bowtie to get the mapping summary
            logger.info(skipCount + " lines were skipped because the read was too short");
            logger.info("</ERROR>");

            int gExitVal = procGenMap.waitFor();
            logger.info("Process exitValue: " + gExitVal);

            brGStdin.close();
            brGStdErr.close();
        } catch (IOException | InterruptedException ex) {
            logger.error("error Bowtie Mapping genome reads\n");
            logger.error(cmdBowtieMapGenomeReads);
            logger.error(ex.toString());
            throw new IOException(STEP_ID_STRING + ": \"error Bowtie Mapping genome reads " + cmdBowtieMapGenomeReads);
        }

        
    }
    
    
    
    
    /**
     * @throws IOException
     * @throws NullPointerException 
     */
    @Override
    public void verifyInputData()  throws IOException, NullPointerException{

        logger.info(STEP_ID_STRING + ": verify input data");        
        this.setPaths();
                
        if(new File(this.getMappingSoftware()).exists() == false){
            logger.error("mapping software not found at location < " + this.getMappingSoftware() +">");
            throw new IOException("mapping software not found at location < " + this.getMappingSoftware() +">");
        }

        //wired. here is to verify the file end with 1.ebwt, but when call bowtie2 abundant just use s directory
//        String pathToAbunBowtieIndex = this.cleanPath(stepInputData.getDataLocations().getGenomeRootFolder()
//                + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_ABUN_DATA_PATH + ".1.ebwt");
//        if(new File(pathToAbunBowtieIndex).exists() == false){
//            logger.error("abundant sequence bowtie index < " + pathToAbunBowtieIndex +"> not found");
//            throw new IOException("abundant sequence bowtie index  < " + pathToAbunBowtieIndex +"> not found");
//        }
        String pathToAbunBowtieIndex = this.cleanPath(getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_ABUN_BOWTIE2_DATA_PATH );
        if(new File(pathToAbunBowtieIndex).isDirectory()== false){
            logger.error("abundant sequence bowtie index < " + pathToAbunBowtieIndex +"> not found");
            throw new IOException("abundant sequence bowtie index  < " + pathToAbunBowtieIndex +"> not found");
        }
        
        String pathToBowtieGenomeIndex = this.cleanPath(getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_BOWTIE2_PATH  + ".1.bt2");
        if(new File(pathToBowtieGenomeIndex).exists() == false){
            logger.error("genome bowtie index < " + pathToBowtieGenomeIndex +"> not found");
            throw new IOException("genome bowtie index  < " + pathToBowtieGenomeIndex +"> not found");
        }
        
                    
        // check the data files
        String fastqFile1;
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            
            //Fastq 1
            if (sampleData.getFastqFile1()==null) {
                logger.error("no Fastq1 file specified");
                throw new IOException("no Fastq1 file specified");
            }
            fastqFile1 = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
                        
            
            if (new File(this.cleanPath(fastqFile1)).exists()==false){
                logger.error("unzipFastqFiles: fastq File1 <" 
                  + fastqFile1 + "> does not exist");
                throw new IOException("unzipFastqFiles: fastq File1 <" 
                  + fastqFile1 + "> does not exist");
            }
            if (fastqFile1.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.error("unzipFastqFiles: incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException("unzipFastqFiles: incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
            
            // we dont check for Fastq 2 as this is single mapping
            
        }
        logger.info("passed");        
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

        configData.put(ID_SOFTWARE, "/usr/bin/bowtie2");
        configData.put(ID_REF_GENOME, "hsa");
        configData.put(ID_MISMATCHES, 1);
        configData.put(ID_ALIGN_MODE, "NODEF");
        configData.put(ID_THREADS, 8);

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
     * @return the AlignMode
     */
    public String getAlignMode() {
        return AlignMode;
    }

    /**
     * @param AlignMode the AlignMode to set
     */
    public void setAlignMode(String AlignMode) {
        this.AlignMode = AlignMode;
    }

    /**
     * @return the NoOfMismatches
     */
    public int getNoOfMismatches() {
        return NoOfMismatches;
    }

    /**
     * @param NoOfMismatches the NoOfMismatches to set
     */
    public void setNoOfMismatches(int NoOfMismatches) {
        this.NoOfMismatches = NoOfMismatches;
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


}

