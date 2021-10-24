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
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;

/**
 * Map reads to contaminants, RNA, Reference Genome and output in SAM format
 *
 * Input is a un-collapsed FASTA file. For mRNA we can't make the same assumptions as smallRNA
 *
 * @author sr
 */
public class StepBowtieMapPairedReads extends NGSStep{

    static Logger logger = LogManager.getLogger();

    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.DATABLE;

    public  static final String STEP_ID_STRING                  = "BowtieMapPairedReads";

    private static final String ID_SOFTWARE                     = "adapter_software";
    private static final String ID_REF_GENOME                   = "host";
    private static final String ID_MISMATCHES                   = "no_of_mismatches";
    private static final String ID_ALIGN_MODE                   = "alignment_mode";
    private static final String ID_THREADS                      = "no_of_threads";

    private static final String INFILE_EXTENSION                = ".trim.clp.fasta";
    private static final String FASTQ_ABUNALN_EXTENSION         = ".trim.clp.abun.fasta";
    private static final String FASTQ_ABUNUNALN_EXTENSION       = ".trim.clp.notabun.fasta";
    private static final String SAM_ABUNALN_EXTENSION           = ".trim.clp.abun.sam";
    private static final String FASTQ_GENALN_EXTENSION          = ".trim.clp.gen.fasta";
    private static final String FASTQ_UNALN_EXTENSION           = ".trim.clp.unmap.fasta";
    private static final String SAM_GENALN_EXTENSION            = ".trim.clp.gen.sam";
    private static final String MAPPING_SUMMARY_EXTENSION       = ".trim.clp.gen.mapping.txt";

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
    
    private  ArrayList<String>  mapGenStdErr;

    public StepBowtieMapPairedReads() {
        classSubtype = NGSStepSubclass.DATABLE;
    }
    
    /**
     *
     * @param sid StepInputData
     *
     */
    public StepBowtieMapPairedReads(InputDataForStep sid) {
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    
    
    @Override
    public String shortStepDescription(){
      return "Performs bowtie2 mapping of paired reads to reference genome";
    }
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "Performs bowtie2 mapping of paired reads to refence genome.\n"
              + "The step requires the Bowtie2 software to be installed.\n"
              + "The location can be specified in the YAML configuration file\n"
              + "along with mapping parameters and reference genome";
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
            throw new NullPointerException("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_REF_GENOME)==null) {
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MISMATCHES)==null) {
            throw new NullPointerException("<" + ID_MISMATCHES + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_ALIGN_MODE)==null) {
            throw new NullPointerException("<" + configData.get(ID_ALIGN_MODE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_REF_GENOME)==null) {
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
        
        
        
        try{
            Integer.parseInt((String) configData.get(ID_THREADS));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_THREADS + " <" + configData.get(ID_THREADS) + "> is not an integer");
        }        
        if (Integer.parseInt((String) configData.get(ID_THREADS)) <= 0){
            throw new IllegalArgumentException(ID_THREADS + " <" + (String) configData.get(ID_THREADS) + "> must be positive");
        }
        this.setNoOfThreads(Integer.parseInt((String) configData.get(ID_THREADS)));

        
        try{
            Integer.parseInt((String) configData.get(ID_MISMATCHES));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_MISMATCHES + " <" + ID_MISMATCHES + "> is not an integer");
        }        
        if (Integer.parseInt((String) configData.get(ID_MISMATCHES)) <= 0){
            throw new IllegalArgumentException(ID_MISMATCHES + " <" + (String) configData.get(ID_MISMATCHES) + "> must be positive");
        }
        this.setNoOfMismatches(Integer.parseInt((String) configData.get(ID_MISMATCHES)));
        

        this.setReferenceGenome((String) configData.get(ID_REF_GENOME));
        if(this.getReferenceGenome().length() !=3 ){
            throw new IllegalArgumentException(ID_REF_GENOME + " <" + (String) configData.get(ID_MISMATCHES) + "> must be a 3 letter string");            
        }
        this.setAlignMode((String) configData.get(ID_ALIGN_MODE));
        this.setMappingSoftware((String) configData.get(ID_SOFTWARE));
        

        logger.info("passed");
    }
    
    
    
    
    
    
    @Override
    public void execute()  throws IOException{

        this.setPaths();
        this.verifyInputData();

        Boolean fA = new File(outFolder).mkdir();
        if (fA) {
            logger.info("created output folder <" + outFolder + "> for results");
        }
        String mappingCmd = this.getMappingSoftware();
        logger.info("Mapping software is " + mappingCmd);
        String pathToBowtieGenomeIndex = getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + ReferenceDataLocations.ID_REL_BOWTIE_PATH;
        
        
        
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()) {
            ArrayList<String> cmd = new ArrayList<>();
            try {
                SampleDataEntry sampleData = (SampleDataEntry) itSD.next();

                /*
                bowtie -a -v 2 e_coli --suppress 1,5,6,7 -c ATGCATCATGCGCCA

                use:
                    -f input files are FASTA format

                    -v option (which ignores quality values) since we are using FASTA files,
                    --best (to order the matches)
                    -m 2 because we only want reads that map to a unique location 
                        (we allow 2, because some miRNAs have two locations)
                    -p 4 (no of threads)
                    --al aligned reads
                    --un unaligned reads
                    --sam SAM file name

                 */
                this.mapAbundantReads(sampleData);
                this.mapReadsToGenome(sampleData);

                /*
                    write out mapping summary
                */
                String mappingOutputFile = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", MAPPING_SUMMARY_EXTENSION);
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
            }
        }

    }

    
    
    /**
     * Maps input reads to the supplied reference abundant sequences
     * 
     * @param sampleData
     * @throws IOException
     * @throws InterruptedException 
     */
    private void mapAbundantReads(SampleDataEntry sampleData) throws IOException, InterruptedException{
        ArrayList<String> cmd = new ArrayList<>();

        String mappingCmd = this.getMappingSoftware();

        String pathToBowtieGenomeIndex = getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + ReferenceDataLocations.ID_REL_BOWTIE_PATH;
        
            cmd.add(mappingCmd);
            String pathToBowtieIndex = this.getRootDataFolder()
                    + FILESEPARATOR + this.getReferenceGenome() + ReferenceDataLocations.ID_REL_ABUN_DATA_PATH;
            cmd.add(pathToBowtieIndex);

            String fastqTrimmedInputFile = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION);
            fastqTrimmedInputFile = fastqTrimmedInputFile.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR).trim();
            cmd.add("-f");
            cmd.add(fastqTrimmedInputFile);

            cmd.add("-v " + this.getAlignMode());
            cmd.add("--best");
            cmd.add("-m 2");

            String fastqAbundantAln = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", FASTQ_ABUNALN_EXTENSION);
            String fastqAbundantUnAln = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", FASTQ_ABUNUNALN_EXTENSION);
            String samAbundantAln = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", SAM_ABUNALN_EXTENSION);
            cmd.add("--al " + fastqAbundantAln);
            cmd.add("--un " + fastqAbundantUnAln);
            cmd.add("--sam " + samAbundantAln);
            cmd.add("-p " + this.getNoOfThreads());
            cmd.add("");

            String cmdBowtieMapAbunReads = StringUtils.join(cmd, " ");
            cmdBowtieMapAbunReads = cmdBowtieMapAbunReads.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR);
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
            ArrayList<String> mapAbunStdErr = new ArrayList<>();
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
        
    }
    
    
    /**
     * map reads that didnt map to Abundant query sequences to the specified reference genome
     * 
     * @param sampleData
     * @throws IOException
     * @throws InterruptedException 
     */
    private void mapReadsToGenome(SampleDataEntry sampleData) throws IOException, InterruptedException{
                 
        String pathToBowtieGenomeIndex = getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + ReferenceDataLocations.ID_REL_BOWTIE_PATH;
                
        ArrayList cmd = new ArrayList<>();
        cmd.add(this.getMappingSoftware());
        cmd.add(pathToBowtieGenomeIndex);

        cmd.add("-f");
        cmd.add(fastqAbundantUnAln);

        cmd.add("-v" + this.getNoOfMismatches());
        cmd.add("--best");
        cmd.add("-m 2");

        fastqGenomeAln = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", FASTQ_GENALN_EXTENSION);
        String fastqGenomeUnAln = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", FASTQ_UNALN_EXTENSION);
        String samGenomeAln = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", SAM_GENALN_EXTENSION);
        cmd.add("--al " + fastqGenomeAln);
        cmd.add("--un " + fastqGenomeUnAln);
        cmd.add("--sam " + samGenomeAln);
        cmd.add("-p " + this.getNoOfThreads());

        String cmdBowtieMapGenomeReads = StringUtils.join(cmd, " ");
        cmdBowtieMapGenomeReads = cmdBowtieMapGenomeReads.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR);
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
        ArrayList<String> mapGenStdErr = new ArrayList<>();
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

        
    }
    
    @Override
    public void verifyInputData()  throws IOException, NullPointerException{
        logger.info("verify input data");
        
        Validate.notNull((String) this.getMappingSoftware());
        
        if (this.getNoOfThreads() <= 0)
        {
            logger.error("number of threads <" + this.getNoOfThreads() + "> must be positive");    
            return;            
        }
                    
        // check the data files
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            String fastqFile1 = (String)sampleData.getFastqFile1();
            String fastqFile2 = (String)sampleData.getFastqFile2();
            
            //Fastq 1
            if (fastqFile1==null) throw new IOException("no Fastq1 file specified");
            
            if ((new File(fastqFile1)).exists()==false){
                throw new IOException("unzipFastqFiles: fastq File1 <" 
                  + fastqFile1 + "> does not exist");
            }
            if (fastqFile1.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                throw new IOException("unzipFastqFiles: incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
            
            
            // dont check for Fastq 2 as this is single mapping
                        
            
        }
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

        configData.put(ID_SOFTWARE, "/usr/local/bin/bowtie");
        configData.put(ID_REF_GENOME, "hsa");
        configData.put(ID_MISMATCHES, 2);
        configData.put(ID_ALIGN_MODE, "v");
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
