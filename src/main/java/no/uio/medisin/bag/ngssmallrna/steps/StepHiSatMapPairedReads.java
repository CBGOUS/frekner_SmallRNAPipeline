/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
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
 * map paired end reads using HiSat package
 * as we are using longer reads compared to smallRNA seq, we skip
 * Input is a collapsed FASTA file
 * 
 * modified from class for bowtie mapping written by sr
 * 
 * @author zxf
 */


public class StepHiSatMapPairedReads extends NGSStep{

    static Logger               logger                          = LogManager.getLogger();

    public  static final NGSStepSubclass STEP_SUBCLASS          = NGSStepSubclass.DATABLE;

    public  static final String CRLF                            = System.getProperty("line.separator");
    public  static final String STEP_ID_STRING                  = "HisatMapPairedReads";

    private static final String ID_SOFTWARE                     = "mappingSoftware";
    private static final String ID_REF_GENOME                   = "host";
    private static final String ID_MISMATCHES                   = "noOfMismatches";
    private static final String ID_ALIGN_MODE                   = "alignmentMode";
    private static final String ID_THREADS                      = "noOfThreads";

    private static final String INFILE_EXTENSION                = ".trim.fastq";
    private static final String FASTQ_ABUNALN_EXTENSION         = ".trim.abun.fastq";
    private static final String FASTQ_ABUNUNALN_EXTENSION       = ".trim.notabun.fastq";
    private static final String SAM_ABUNALN_EXTENSION           = ".trim.abun.sam";
    private static final String FASTQ_PAIR_ALN_EXTENSION        = ".trim.gen.pair.aln.fastq";
    private static final String FASTQ_PAIR_UNCONC_EXTENSION     = ".trim.pairunc.unconc.fastq";
    private static final String FASTQ_UNPAIR_UNALN_EXTENSION    = ".trim.unpair.unaln.fastq";
    private static final String FASTQ_UNPAIR_UNCONC_EXTENSION   = ".trim.unpair.unconc.fastq";
    private static final String SAM_GENALN_EXTENSION            = ".trim.clp.gen.sam";
    private static final String MAPPING_SUMMARY_EXTENSION       = ".trim.clp.gen.mapping.txt";
    
    private static final String RAW_INPUT_EXTENSION             = ".fastq.gz";

    private             String  mappingSoftware                 = "";
    private             String  AlignMode                       = "";
    private             int     NoOfMismatches                  = 2;
    private             int     NoOfThreads                     = 4;
    private             String  rootDataFolder                  = "";
    private             String  ReferenceGenome                 = "";
    
    private             String  fastqTrimmedInputFile1           = "";
    private             String  fastqTrimmedInputFile2           = "";
    private             String  fastqAbundantAln                = "";
    private             String  fastqAbundantUnAln              = "";
    private             String  fastqGenomeUnpairUnAln          = "";
    private             String  fastqGenomeUnpairUnConc         = "";
    private             String  fastqGenomePairUnConc           = "";
    private             String  fastqGenomePairAlConc           = "";
    private             String  fastqGenomeUnAln                = "";
    private             String  samGenomeAln                    ="";
    
    private  ArrayList<String>  mapGenStdErr;
    private  ArrayList<String>  mapAbunStdErr;
    
    
    public StepHiSatMapPairedReads() {
        classSubtype = NGSStepSubclass.DATABLE;
    }

    /**
     *
     * @param sid StepInputData
     *
     */
    public StepHiSatMapPairedReads(InputDataForStep sid) {
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
   
  
    @Override
    public String shortStepDescription(){
      return "map paired end reads using HiSat package";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "map paired end reads using HiSat package.\n"
              + "The step requires the HiSat software to be installed.\n"
              + "The location can be specified in the YAML configuration file.\n"
              + "Input is a collapsed FASTA file";
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
        if(configData.get(ID_THREADS)==null) {
            logger.error("<" + ID_THREADS + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_THREADS + "> : Missing Definition in Configuration File");
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
        String hisatCmdFile = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + ".hisat.sh";
        logger.info("hisat command file is <" + hisatCmdFile + ">");
        BufferedWriter bwHC = new BufferedWriter(new FileWriter(new File(hisatCmdFile)));
        
        
        String mappingCmd = this.getMappingSoftware();
        logger.info("Mapping software is " + mappingCmd);
        
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()) {
            ArrayList<String> cmd = new ArrayList<>();
            try {
                SampleDataEntry sampleData = (SampleDataEntry) itSD.next();

                bwHC.write(this.mapReadsToGenome(sampleData) + "\n");

                /*
                    write out mapping summary
                String mappingOutputFile = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, MAPPING_SUMMARY_EXTENSION);
                BufferedWriter bwMO = new BufferedWriter(new FileWriter(new File(mappingOutputFile)));
                for (String mapLine : mapGenStdErr) {
                    bwMO.write(mapLine + "\n");
                }
                bwMO.write("\n\n" + "+" + StringUtils.repeat("-", 60) + "+" + "\n");
                bwMO.write("original FASTQ1 source" + sampleData.getFastqFile1() + "\n");
                bwMO.write("original FASTQ2 source" + sampleData.getFastqFile2() + "\n");
                bwMO.write(fastqTrimmedInputFile1 + "\n");

                bwMO.write(fastqGenomeUnpairUnAln + "\n");
                bwMO.write(fastqGenomeUnpairUnConc + "\n");
                bwMO.write(fastqGenomePairUnConc + "\n");
                bwMO.write(fastqGenomePairAlConc + "\n");

                bwMO.write("\n\n" + "+" + StringUtils.repeat("-", 60) + "+" + "\n");

                bwMO.close();
                */
            } catch (IOException | InterruptedException ex) {
                logger.error("error executing Bowtie Mapping command\n");
                logger.error(cmd);
                logger.error(ex.toString());
                throw new IOException(STEP_ID_STRING + ": \"error executing Bowtie Mapping command " + cmd);
            }
        }
        logger.info(STEP_ID_STRING + ": completed");
        bwHC.close();
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

            fastqTrimmedInputFile1 = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
            cmd.add("-f");
            cmd.add(fastqTrimmedInputFile1);

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
    private String mapReadsToGenome(SampleDataEntry sampleData) throws IOException, InterruptedException{
                 
        logger.info(STEP_ID_STRING + ": mapping reads to genome");
        String cmdHiSatMapGenomeReads = "";
        //try{
            String pathToHisatIndex = this.cleanPath(getStepInputData().getDataLocations().getGenomeRootFolder()
                    + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_HISAT_PATH);

            ArrayList cmd = new ArrayList<>();
            cmd.add(this.getMappingSoftware());
            cmd.add("-x");
            cmd.add(pathToHisatIndex);

            cmd.add("-q");
            fastqTrimmedInputFile1 = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
            cmd.add("-1 " + fastqTrimmedInputFile1);
            fastqTrimmedInputFile2 = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile2().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
            cmd.add("-2 " + fastqTrimmedInputFile2);

//            cmd.add("-" + this.getAlignMode());
//            cmd.add("--best");
            
            cmd.add("-p " + this.getNoOfThreads());

            fastqGenomeUnpairUnConc = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, FASTQ_UNPAIR_UNCONC_EXTENSION));
            fastqGenomeUnpairUnAln = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, FASTQ_UNPAIR_UNALN_EXTENSION));
            fastqGenomePairAlConc = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, FASTQ_PAIR_ALN_EXTENSION));
            fastqGenomePairUnConc = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, FASTQ_PAIR_UNCONC_EXTENSION));

            cmd.add("--al " + fastqGenomeUnpairUnConc);
            cmd.add("--un " + fastqGenomeUnpairUnAln);
            cmd.add("--un-conc " + fastqGenomePairUnConc);
            cmd.add("--al-conc " + fastqGenomePairAlConc);
            
            
            samGenomeAln = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, SAM_GENALN_EXTENSION));
            cmd.add("-S " + samGenomeAln);
            

            cmdHiSatMapGenomeReads = this.cleanPath(StringUtils.join(cmd, " "));
            logger.info("HiSat Map Genome Reads command:\t" + cmdHiSatMapGenomeReads);
            ArrayList cmdSet = new ArrayList<>();
            cmdSet.add("#  " + Paths.get(fastqTrimmedInputFile1).getFileName().toString() 
                    + " | " + (Paths.get(fastqTrimmedInputFile2)).getFileName().toString());
            cmdSet.add("pigz -d -p 4 " + fastqTrimmedInputFile1);
            cmdSet.add("pigz -d -p 4 " + fastqTrimmedInputFile2);
            cmdSet.add(cmdHiSatMapGenomeReads);
            cmdSet.add("pigz -p 4 " + fastqTrimmedInputFile1);
            cmdSet.add("pigz -p 4 " + fastqTrimmedInputFile2);
            cmdSet.add("pigz -p 4 " + fastqGenomeUnpairUnConc);
            cmdSet.add("pigz -p 4 " + fastqGenomeUnpairUnAln);
            cmdSet.add("pigz -p 4 " + fastqGenomePairAlConc);
            cmdSet.add("pigz -p 4 " + fastqGenomePairUnConc);
            cmdSet.add("pigz -p 4 " + samGenomeAln);
            cmdSet.add("\n");
            
            return StringUtils.join(cmdSet, CRLF);
            /*
            Runtime rtGenMap = Runtime.getRuntime();
            Process procGenMap = rtGenMap.exec(cmdHiSatMapGenomeReads);
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
            logger.error(cmdHiSatMapGenomeReads);
            logger.error(ex.toString());
            throw new IOException(STEP_ID_STRING + ": \"error Bowtie Mapping genome reads " + cmdHiSatMapGenomeReads);
        }
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
                
        
        if(new File(this.getMappingSoftware()).exists() == false){
            logger.error("mapping software not found at location < " + this.getMappingSoftware() +">");
            throw new IOException("mapping software not found at location < " + this.getMappingSoftware() +">");
        }


        /*
            we don't map to abundant for now since we can filter these out later using a GTF file 
        String pathToAbunBowtieIndex = this.cleanPath(getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_ABUN_BOWTIE2_DATA_PATH );
        if(new File(pathToAbunBowtieIndex).isDirectory()== false){
            logger.error("abundant sequence bowtie index < " + pathToAbunBowtieIndex +"> not found");
            throw new IOException("abundant sequence bowtie index  < " + pathToAbunBowtieIndex +"> not found");
        }
        */
        String pathToHisatGenomeIndexShortInt =  this.cleanPath(getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_HISAT_PATH  + ".1.bt2");
        
        String pathToHisatGenomeIndexLongInt =  this.cleanPath(getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_HISAT_PATH  + ".1.bt2l");
        if(new File(pathToHisatGenomeIndexShortInt).exists() == false && new File(pathToHisatGenomeIndexLongInt).exists() == false){
            logger.error("no genome hisat short or long index < " 
                    + pathToHisatGenomeIndexShortInt 
                    + "||" 
                    + pathToHisatGenomeIndexLongInt +"> was found");
            throw new IOException("no genome hisat short or long index  < " 
                    + pathToHisatGenomeIndexShortInt 
                    + "||" + pathToHisatGenomeIndexLongInt 
                    +"> was found");
        }        
        
        


                    /*
        // check the data files
        String fastqFile1;
        String fastqFile2;
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
                logger.error("HisatMapPairedReads: fastq File1 <" 
                  + fastqFile1 + "> does not exist");
                throw new IOException("HisatMapPairedReads: fastq File1 <" 
                  + fastqFile1 + "> does not exist");
            }
            if (fastqFile1.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.error("HisatMapPairedReads: incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException("HisatMapPairedReads: incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
            
            //Fastq2
            if (sampleData.getFastqFile2()==null) {
                logger.error("no Fastq1 file specified");
                throw new IOException("no Fastq1 file specified");
            }
            fastqFile2 = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile2().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
                        
            
            if (new File(this.cleanPath(fastqFile2)).exists()==false){
                logger.error("HisatMapPairedReads: fastq File2 <" 
                  + fastqFile2 + "> does not exist");
                throw new IOException("HisatMapPairedReads: fastq File2 <" 
                  + fastqFile2 + "> does not exist");
            }
            if (fastqFile2.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.error("HisatMapPairedReads: incorrect file extension for input file <" 
                  + fastqFile2 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException("HisatMapPairedReads: incorrect file extension for input file <" 
                  + fastqFile2 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
            
        }
        */
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

        configData.put(ID_SOFTWARE, "/usr/local/bin/hisat");
        configData.put(ID_REF_GENOME, "h38");
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

