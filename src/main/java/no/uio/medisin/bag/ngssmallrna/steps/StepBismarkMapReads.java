/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;




/**
 *  Perform mapping of single end BS reads using Bismark (a 3 letter aligner)
 *  the calculation of methylation ratios is performed in a separate step
 * 
 *   Input is a FASTQ file, raw or trimmed
 *   Output is a alignment file in SAM format
 * 
 * 
 * @author sr
 */

public class StepBismarkMapReads extends NGSStep{
    
    private static Logger               logger                      = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS               = NGSStepSubclass.DATABLE;

    public static final String          STEP_ID_STRING              = "BismarkMapSingleReads";

    private static final String         ID_SOFTWARE                 = "pathToBismark";
    private static final String         ID_PATH_TO_BOWTIE           = "pathToBowtie";
    private static final String         ID_PATH_TO_SAMTOOLS         = "pathToSAMTools";
    private static final String         ID_HOST                     = "host";

    private static final String         ID_SKIP_FIRST_NREADS        = "skipFirstNReads";
    private static final String         ID_ALIGN_FIRST_NREADS       = "alignFirstNReads";
    private static final String         ID_MAX_NUM_SEED_MISMATCHES  = "maxNumMismatches";
    private static final String         ID_MAX_SEED_LENGTH          = "maxSeedLength";
    private static final String         ID_MAX_QUALVAL_TOTAL        = "maxQualValTotal";
    private static final String         ID_CHUNK_MEM_SIZE           = "chunkMemoryInMB";
    private static final String         ID_THREADS                  = "noOfThreads";

    private static final String         ID_WRITE_BISMARK_OUTPUT     = "writeBismarkOutput";
    private static final String         ID_WRITE_UNALIGNED_READS    = "writeUnalignedReads";
    private static final String         ID_WRITE_AMBIGUOUS_READS    = "writeAmbiguousReads";
    private static final String         ID_WRITE_SAM_FORMAT         = "writeSAMFormat";

    
    
    private static final String         RAW_INPUT_EXTENSION         = ".fastq.gz";
    private static final String         INFILE_EXTENSION            = ".trim.fastq";
    private static final String         OUTFILE_EXTENSION           = ".trim.bsmap.sam";

    

    private             String          pathToBismark               = "";
    private             String          pathToBowtie                = "";
    private             String          pathToSAMTools               = "";
    private             String          ReferenceGenome             = "";

    private             int             skipFirstNReads             = 1;
    private             int             alignFirstNReads            = -1;
    private             int             maxSeedLength               = 16;
    private             int          maxNumMismatches            = 2;
    private             int             maxQualValTotal             = 70;
    private             int             chunkMemoryInMB             = 0;
    private             int             noOfThreads                 = 4;

    private             Boolean         writeUnalignedReads         = false;
    private             Boolean         writeAmbiguousReads         = false;
    private             Boolean         writeSAMFormat              = false;
    private             Boolean         writeBismarkOutput          = false;

     

    public StepBismarkMapReads(){
        classSubtype = NGSStepSubclass.DATABLE;
    }
    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepBismarkMapReads(InputDataForStep sid){
        classSubtype = NGSStepSubclass.DATABLE;
       stepInputData = sid;
    }
    
    

    @Override
    public String shortStepDescription(){
      return "Perform mapping of single end BS reads using Bismark (a 3 letter aligner)";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "Perform mapping of single end BS reads using Bismark (a 3 letter aligner)\n"
              + "the calculation of methylation ratios is performed in a separate step.\n\n"
              + "Input is a FASTQ file, raw or trimmed\n" 
              + "Output is a alignment file in SAM format";

    }


    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    
    /**
     * This parses out the hashmap containing the run parameters for this step
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{

        logger.info(STEP_ID_STRING + ": verify configuration data");
        
        /*
            the following parameters have to be specfied in the configuration file for the 
            analysis to proceed
        */
        if(configData == null || configData.get(ID_SOFTWARE)==null) {
            logger.info("<" + ID_SOFTWARE + "> : Missing Definition in Configuration File");
            logger.error("<" + ID_SOFTWARE + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_SOFTWARE + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_HOST)==null) {
            logger.info("<" + configData.get(ID_HOST) + "> : Missing Definition in Configuration File");
            logger.error("<" + configData.get(ID_HOST) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_HOST) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_THREADS)==null) {
            logger.info("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
            logger.error("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
        }

        
        this.setReferenceGenome((String) configData.get(ID_HOST));
        if(this.getReferenceGenome().length() !=3 ){
            logger.error(ID_HOST + " <" + configData.get(ID_HOST) + "> must be a 3 letter string");            
            throw new IllegalArgumentException(ID_HOST + " <" + configData.get(ID_HOST) + "> must be a 3 letter string");            
        }
        this.setPathToBismarkMap((String) configData.get(ID_SOFTWARE));
        this.setPathToBowtie((String) configData.get(ID_PATH_TO_BOWTIE));
        this.setPathToSAMTools((String) configData.get(ID_PATH_TO_SAMTOOLS));
        

        String chk;

        chk = checkParameter("Integer", ID_MAX_NUM_SEED_MISMATCHES, Integer.toString((Integer)configData.get(ID_MAX_NUM_SEED_MISMATCHES)), "0", "1", logger);
        if(chk!=null)
            this.setMaxNumMismatches((Integer)configData.get(ID_MAX_NUM_SEED_MISMATCHES));
            
        chk = checkParameter("Integer", ID_MAX_SEED_LENGTH, Integer.toString((Integer)configData.get(ID_MAX_SEED_LENGTH)), "0", "NA", logger);
        if(chk!=null)
            this.setMaxSeedLength((Integer)configData.get(ID_MAX_SEED_LENGTH));
                        
        chk = checkParameter("Integer", ID_SKIP_FIRST_NREADS, Integer.toString((Integer)configData.get(ID_SKIP_FIRST_NREADS)), "0", "NA", logger);
        if(chk!=null)
            this.setSkipFirstNReads((Integer)configData.get(ID_SKIP_FIRST_NREADS));
            
        chk = checkParameter("Integer", ID_ALIGN_FIRST_NREADS, Integer.toString((Integer)configData.get(ID_ALIGN_FIRST_NREADS)), "0", "NA", logger);
        if(chk!=null)
            this.setAlignFirstNReads((Integer)configData.get(ID_ALIGN_FIRST_NREADS));
                        
        chk = checkParameter("Integer", ID_THREADS, Integer.toString((Integer)configData.get(ID_THREADS)), "0", "12", logger);
        if(chk!=null)
            this.setNoOfThreads((Integer)configData.get(ID_THREADS));
            
        chk = checkParameter("Integer", ID_MAX_QUALVAL_TOTAL, Integer.toString((Integer)configData.get(ID_MAX_QUALVAL_TOTAL)), "0", "NA", logger);
        if(chk!=null)
            this.setMaxQualValTotal((Integer)configData.get(ID_MAX_QUALVAL_TOTAL));
            
        chk = checkParameter("Integer", ID_CHUNK_MEM_SIZE, Integer.toString((Integer)configData.get(ID_CHUNK_MEM_SIZE)), "512", "NA", logger);
        if(chk!=null)
            this.setChunkMemoryInMB((Integer)configData.get(ID_CHUNK_MEM_SIZE));
            

        chk = checkParameter("Boolean", ID_WRITE_BISMARK_OUTPUT, Boolean.toString((Boolean)configData.get(ID_WRITE_BISMARK_OUTPUT)), "NA", "NA", logger);
        if(chk!=null){
            this.setWriteBismarkOutput((Boolean)configData.get(ID_WRITE_BISMARK_OUTPUT));
        }
            
        chk = checkParameter("Boolean", ID_WRITE_SAM_FORMAT, Boolean.toString((Boolean)configData.get(ID_WRITE_SAM_FORMAT)), "NA", "NA", logger);
        if(chk!=null){
            this.setWriteSAMFormat((Boolean)configData.get(ID_WRITE_SAM_FORMAT));
        }
            
        chk = checkParameter("Boolean", ID_WRITE_UNALIGNED_READS, Boolean.toString((Boolean)configData.get(ID_WRITE_UNALIGNED_READS)), "NA", "NA", logger);
        if(chk!=null){
            this.setWriteUnalignedReads((Boolean)configData.get(ID_WRITE_UNALIGNED_READS));
        }
            
        chk = checkParameter("Boolean", ID_WRITE_AMBIGUOUS_READS, Boolean.toString((Boolean)configData.get(ID_WRITE_AMBIGUOUS_READS)), "NA", "NA", logger);
        if(chk!=null){
            this.setWriteAmbiguousReads((Boolean)configData.get(ID_WRITE_AMBIGUOUS_READS));
        }
            

        logger.info("passed");
    }
    
    
    
    /**
     * bismark map single end reads
     * 
     * 
     * @throws IOException 
     */
    @Override
    public void execute() throws IOException{
        /* example call.:
        -C or --cytosine contexts acquired
        -sm or --sequencing mode
        -nt or --num threads
        -I or --input file        
        -D or --dbsnp (dbsnp VCF file)
        -R or --reference sequence 
        -L or --intervals < intervals file name >
        -vfn1 or --vcf file name 1 CpG sites
        -vfn2 or --vcf file name 2 SNP sites
        -out modes or --output modes (controls output of sites using threshold)
        
        -cpgreads (not sure I completely understand this parameter)
        -stand call conf < standard min confidence threshold for calling >
        -stand emit conf < standard min confidence threshold for emitting >
        -mmq,--min mapping quality score < min mapping quality score >
        -mbq,--min base quality score < min base quality score >
        -minConv,--minmum cytosine converted < minmum number of cytosine converted >
        4.1
        - Add read group tag to BAM file
        - Indel realignment
            Find indel region
            Realign in the indel region
        - Mark duplicated reads
        - Base quality recalibration
            Count Covariant
            Write recalibrated base quality score into BAM file
            Re-Count Covariant
            Generate recalibration plot
        - Bis-SNP genotyping
        - Filter fake SNPs
        - Generate bed file or wig file for SNP/DNA methylation visualization
            Converted to .wig format
            Converted to .bed format
            Converted to .bedgraph format
            Extract cytosine coverage information to .bedgraph format
        
        GATK engine requires BAM file to have ReadGroup tag. 
        When your own BAM file does not contained Read group tag. 
        Download Picard tools http://sourceforge.net/projects/picard/files/, t
        hen using the following command to add Read group tag to BAM file.
        
            java -Xmx4g -jar AddOrReplaceReadGroups.jar I=sample.withoutRG.bam O=sample.withRG.bam ID=readGroup_name
            LB=readGroup_name PL=illumina PU=run SM=sample_name CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
            SORT_ORDER=coordinate
        
            java -Xmx4g -jar 
            BisSNP-0.71.jar 
            -T BisulfiteGenotyper 
            -R hg18_unmasked.plusContam.fa 
            -D dbsnp_135.hg18.sort.vcf 
            -I normalMerge_chr11-7M-9M.nodups.withhead.bam 
            -vfn1 cpg.raw.vcf -vfn2 snp.raw.vcf 
            -stand_call_conf 30 
            -stand_emit_conf 0 
            -L Interval_list-test.bed        
        */
        logger.info(STEP_ID_STRING + ": execute step");        
        
        Boolean fA = new File(outFolder).mkdir();
        if (fA) {
            logger.info("created output folder <" + outFolder + "> for results");
        }
        
        
        String fastqFile1in = "";
        String fastqFile1out = "";
        String fastqFile2in = "";
        String fastqFile2out = "";
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();                
                fastqFile1in = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION);
                fastqFile1out = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, OUTFILE_EXTENSION);

    
                if (new File(fastqFile1in).exists()==false){
                    logger.error(STEP_ID_STRING + ": fastq File 1 <" + fastqFile1in + "> does not exist");
                    throw new IOException(STEP_ID_STRING + ": fastq File 1 <" + fastqFile1in + "> does not exist");
                }
                    
                String cmdBismarkMap = "";   
                ArrayList<String> cmd = new ArrayList<>();
                cmd.add(this.getPathToBismarkMap() + " -fastq");
                
                cmd.add("--single_end " + fastqFile1in);
                cmd.add("-N " + this.getMaxNumMismatches());
                cmd.add("-L " + this.getMaxSeedLength());
                cmd.add("-chunkmbs " + this.getChunkMemoryInMB());
                cmd.add("-multicore " + this.getNoOfThreads());
                
                
                cmd.add("--samtools_path " + this.getPathToSAMTools());
                cmd.add("--path_to_bowtie " + this.getPathToBowtie());
                cmd.add("--bowtie2");
                cmd.add("-o " + outFolder);
                
                
                if(this.getWriteUnalignedReads())
                    cmd.add("--un");
                



                
                String pathToGenomeFA = this.cleanPath(getStepInputData().getDataLocations().getGenomeRootFolder()
                        + FILESEPARATOR + this.getReferenceGenome() + "/Sequence/WholeGenomeFasta" + FILESEPARATOR);
                cmd.add(pathToGenomeFA);

                    
                                
                cmdBismarkMap = this.cleanPath(StringUtils.join(cmd, " "));
                logger.info("BSMap command is: " + cmdBismarkMap);
                /*
                Runtime rt = Runtime.getRuntime();
                Process proc = rt.exec(cmdBismarkMap);
                ReadStdIOStream inputStream = new ReadStdIOStream("stdin", proc.getInputStream(), STEP_ID_STRING);
                ReadStdIOStream errorStream = new ReadStdIOStream("stderr", proc.getErrorStream(), STEP_ID_STRING);
                inputStream.start();
                errorStream.start();
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
                
                ArrayList<String>  mapAbunStdErr;

                mapAbunStdErr = new ArrayList<>();
                while ((line = brAStdErr.readLine()) != null) {
                    mapAbunStdErr.add(line);
                }
                

                logger.info("</ERROR>");

                int exitVal = proc.waitFor();
                logger.info("Process exitValue: " + exitVal);

                brAStdin.close();
                brAStdErr.close();
                logger.info(STEP_ID_STRING + ": done");

                String samAlnOutput = this.cleanPath(outFolder + FILESEPARATOR 
                    + sampleData.getFastqFile1().split("_")[0].trim()+".bsmap_output");
                try(BufferedWriter bwAO = new BufferedWriter(new FileWriter(new File(samAlnOutput)))){
                    for(String logString:mapAbunStdErr){
                        bwAO.write(logString + "\n");
                    }
                }
                catch(IOException exIO){
                    logger.error("error writing BSMap summary file <" + samAlnOutput +  ">");
                    logger.info("error writing BSMap summary file <" + samAlnOutput +  ">");
                    throw new IOException("error writing BSMap summary file <" + samAlnOutput +  ">" + exIO);
                }
                */                
                
            }
//            catch(IOException | InterruptedException ex ){
            catch(IOException ex ){
                logger.info("error executing bsmap command:\n" + ex.toString());
                logger.error("error executing bsmap command:\n" + ex.toString());
                throw new IOException("error executing bsmap command");
            }
                
        }
        
        logger.info(STEP_ID_STRING + ": completed");
    }
    
    
    
            
    /**
     * this should be called prior to executing the step.
     * check unzip software exists and input files are available
     * 
     * @throws IOException
     */
    @Override
    public void verifyInputData() throws IOException{
        
        logger.info("verify input data");        
        this.setPaths();
        
        if(new File(this.getPathToBismarkMap()).exists() == false){
            logger.info("Bismark software not found at location < " + this.getPathToBismarkMap() +">");
            logger.error("BSMap software not found at location < " + this.getPathToBismarkMap() +">");
            throw new IOException("BSMap software not found at location < " + this.getPathToBismarkMap() +">");
        }
                
        if(new File(this.getPathToBowtie()).exists() == false){
            logger.info("Bowtie software not found at location < " + this.getPathToBowtie() +">");
            logger.error("Bowtie software not found at location < " + this.getPathToBowtie() +">");
            throw new IOException("Bowtie software not found at location < " + this.getPathToBowtie() +">");
        }
                
        if(new File(this.getPathToSAMTools()).exists() == false){
            logger.info("SAMTools software not found at location < " + this.getPathToSAMTools() +">");
            logger.error("SAMTools software not found at location < " + this.getPathToSAMTools() +">");
            throw new IOException("SAMTools software not found at location < " + this.getPathToSAMTools() +">");
        }
                

                
                            
        // check the data files
        String fastqFile1in = "";
        String fastqFile1out = "";
        for(SampleDataEntry sampleData: this.getStepInputData().getSampleData()){
            
            
            //Fastq 1
            if (sampleData.getFastqFile1()==null) {
                logger.error("no Fastq1 file specified");
                throw new IOException("no Fastq1 file specified");
            }
            String fastqFile1 = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION);
            
/*
            if (new File(fastqFile1).exists()==false){
                logger.error(STEP_ID_STRING + ": fastq 1 File <" + fastqFile1 + "> does not exist");
                throw new IOException(STEP_ID_STRING + ": fastq 1 File <" + fastqFile1 + "> does not exist");
            }
            
            
            //Fastq 2
            if (sampleData.getFastqFile2()==null) continue;
            String fastqFile2in = "";
            String fastqFile2out = "";

            
            if ((new File(fastqFile2in)).exists()==false ){
                logger.error(STEP_ID_STRING + ": fastq 2 File <" + fastqFile2in + "> does not exist");
                throw new IOException(STEP_ID_STRING + ": fastq 2 File <" + fastqFile2in + "> does not exist");
            }
*/                        
        }
        logger.info("passed");

    }
    
    
    
    @Override
    public void verifyOutputData(){
        logger.info("no output verification required");
    }
    
    
    
    @Override
    public NGSStepSubclass getStepSubclass(){
        return STEP_SUBCLASS;
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
        
        configData.put(ID_SOFTWARE, "/usr/local/bin/bismark");
        configData.put(ID_PATH_TO_BOWTIE, "/usr/local/bin/bowtie");
        configData.put(ID_PATH_TO_SAMTOOLS, "/usr/local/bin/samtools");
        configData.put(ID_HOST, "hsa");
        
        configData.put(ID_ALIGN_FIRST_NREADS, 0);
        configData.put(ID_SKIP_FIRST_NREADS, 0);
        
        configData.put(ID_MAX_NUM_SEED_MISMATCHES, 2);
        configData.put(ID_MAX_SEED_LENGTH, 28);
        configData.put(ID_THREADS, 4);
        configData.put(ID_MAX_QUALVAL_TOTAL, 70);
        configData.put(ID_CHUNK_MEM_SIZE, 512);
        
        configData.put(ID_WRITE_BISMARK_OUTPUT, false);
        configData.put(ID_WRITE_UNALIGNED_READS, true);
        configData.put(ID_WRITE_AMBIGUOUS_READS, true);
        configData.put(ID_WRITE_SAM_FORMAT, true);


        return configData;
        
    }

    
    
    
    
    

    /**
     * @return the pathToBSMap
     */
    public String getPathToBismarkMap() {
        return pathToBismark;
    }

    /**
     * @param pathToBismarkMap the pathToBSMap to set
     */
    public void setPathToBismarkMap(String pathToBismarkMap) {
        this.pathToBismark = pathToBismarkMap;
    }

    /**
     * @return the noOfMismatches
     */
    public int getMaxNumMismatches() {
        return maxNumMismatches;
    }

    /**
     * @param maxNumMismatches the noOfMismatches to set
     */
    public void setMaxNumMismatches(int maxNumMismatches) {
        this.maxNumMismatches = maxNumMismatches;
    }

    /**
     * @return the noOfThreads
     */
    public int getNoOfThreads() {
        return noOfThreads;
    }

    /**
     * @param noOfThreads the noOfThreads to set
     */
    public void setNoOfThreads(int noOfThreads) {
        this.noOfThreads = noOfThreads;
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
     * @return the qualityThresholdTrimValue
     */
    public int getChunkMemoryInMB() {
        return chunkMemoryInMB;
    }

    /**
     * @param chunkMemoryInMBValue the qualityThresholdTrimValue to set
     */
    public void setChunkMemoryInMB(int chunkMemoryInMBValue) {
        this.chunkMemoryInMB = chunkMemoryInMBValue;
    }

    /**
     * @return the maxInsertSize
     */
    public int getMaxQualValTotal() {
        return maxQualValTotal;
    }

    /**
     * @param maxQualValTotal the maxInsertSize to set
     */
    public void setMaxQualValTotal(int maxQualValTotal) {
        this.maxQualValTotal = maxQualValTotal;
    }

    /**
     * @return the reportUnmappedReads
     */
    public Boolean getWriteUnalignedReads() {
        return writeUnalignedReads;
    }

    /**
     * @param writeUnalignedReads the reportUnmappedReads to set
     */
    public void setWriteUnalignedReads(Boolean writeUnalignedReads) {
        this.writeUnalignedReads = writeUnalignedReads;
    }

    /**
     * @return the startAtThisRead
     */
    public int getSkipFirstNReads() {
        return skipFirstNReads;
    }

    /**
     * @param skipFirstNReads the startAtThisRead to set
     */
    public void setSkipFirstNReads(int skipFirstNReads) {
        this.skipFirstNReads = skipFirstNReads;
    }

    /**
     * @return the endtAtThisRead
     */
    public int getAlignFirstNReads() {
        return alignFirstNReads;
    }

    /**
     * @param alignFirstNReads the endtAtThisRead to set
     */
    public void setAlignFirstNReads(int alignFirstNReads) {
        this.alignFirstNReads = alignFirstNReads;
    }

    /**
     * @return the maxSeedLength
     */
    public int getMaxSeedLength() {
        return maxSeedLength;
    }

    /**
     * @param maxSeedLength the maxSeedLength to set
     */
    public void setMaxSeedLength(int maxSeedLength) {
        this.maxSeedLength = maxSeedLength;
    }

    /**
     * @return the pathToBowtie
     */
    public String getPathToBowtie() {
        return pathToBowtie;
    }

    /**
     * @param pathToBowtie the pathToBowtie to set
     */
    public void setPathToBowtie(String pathToBowtie) {
        this.pathToBowtie = pathToBowtie;
    }

    /**
     * @return the pathToSAMTools
     */
    public String getPathToSAMTools() {
        return pathToSAMTools;
    }

    /**
     * @param pathToSAMTools the pathToSAMTools to set
     */
    public void setPathToSAMTools(String pathToSAMTools) {
        this.pathToSAMTools = pathToSAMTools;
    }

    /**
     * @return the writeAmbiguousReads
     */
    public Boolean getWriteAmbiguousReads() {
        return writeAmbiguousReads;
    }

    /**
     * @param writeAmbiguousReads the writeAmbiguousReads to set
     */
    public void setWriteAmbiguousReads(Boolean writeAmbiguousReads) {
        this.writeAmbiguousReads = writeAmbiguousReads;
    }

    /**
     * @return the writeSAMFormat
     */
    public Boolean getWriteSAMFormat() {
        return writeSAMFormat;
    }

    /**
     * @param writeSAMFormat the writeSAMFormat to set
     */
    public void setWriteSAMFormat(Boolean writeSAMFormat) {
        this.writeSAMFormat = writeSAMFormat;
    }

    /**
     * @return the writeBismarkOutput
     */
    public Boolean getWriteBismarkOutput() {
        return writeBismarkOutput;
    }

    /**
     * @param writeBismarkOutput the writeBismarkOutput to set
     */
    public void setWriteBismarkOutput(Boolean writeBismarkOutput) {
        this.writeBismarkOutput = writeBismarkOutput;
    }



    
    
    
    
    
}
