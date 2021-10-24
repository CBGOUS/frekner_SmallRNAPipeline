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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import no.uio.medisin.bag.core.sequence.SimpleSeq;
import no.uio.medisin.bag.ngssmallrna.pipeline.ReferenceDataLocations;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;




/**
 *  Perform mapping of BS reads using BSMap
 *  the calculation of methylation ratios is performed in a separate step
 * 
 *   Input is a FASTQ file, raw or trimmed
 *   Output is a alignment file 
 * 
 * 
 * @author sr
 */

public class StepBSMapReads extends NGSStep{
    
    private static Logger               logger                      = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS               = NGSStepSubclass.DATABLE;

    public static final String          STEP_ID_STRING              = "BSMapReads";

    private static final String         ID_SOFTWARE                 = "pathToBSmap";
    private static final String         ID_HOST                     = "host";
    private static final String         ID_MISMATCHES               = "noOfMismatches";
    private static final String         ID_SEED_SIZE                = "seedSize";
    private static final String         ID_THREEMAPPING             = "threeNucleoideMapping";
    private static final String         ID_GAPSIZE                  = "gapSize";
    private static final String         ID_KMERCUTOFF               = "kmerCutoffRatio";
    private static final String         ID_REPORT_REPEATS           = "reportRepeatsValue";
    private static final String         ID_QUAL_THRESHOLD_VALUE     = "qualityThresholdTrimValue";
    private static final String         ID_LOW_QUAL_VALUE           = "lowQualityFilterValue";
    private static final String         ID_THREADS                  = "noOfThreads";
    private static final String         ID_MAX_INSERT_SIZE          = "maxInsertSize";
    private static final String         ID_MIN_INSERT_SIZE          = "minInsertSize";
    private static final String         ID_MAP_FIRST_N_NUCS         = "mapFirstNucleotides";
    private static final String         ID_GENOME_INDEX_INTERVAL    = "genomeIndexInterval";
    private static final String         ID_ADAPTERSEQ               = "adapterSequence";
    private static final String         ID_TRIM_ADAPTERSEQ          = "trimAdapterSequence";
    private static final String         ID_INCLUDE_REF_SEQ          = "includeRefSeq";
    private static final String         ID_SKIP_SAM_HEADER          = "skipSAMHeader";
    private static final String         ID_REPORT_UNMAPPED_READS    = "reportUnmappedReads";
    private static final String         ID_START_AT_READ_N          = "startAtThisRead";
    private static final String         ID_END_AT_READ_N            = "endtAtThisRead";
    private static final String         ID_DIGESTIONSITE            = "digestionSite";
    private static final String         ID_RANDOMSEED               = "randomSeed";
    private static final String         ID_MAPPINGSTRAND            = "mappingStrand";
    private static final String         ID_TRANSITIONMAP            = "transitionMap";
    private static final String         ID_MESSAGELEVEL             = "messageLevel";
    
    
    private static final String         INFILE_EXTENSION_TRIM       = ".fastq";
    private static final String         INFILE_EXTENSION_NOTRIM     = ".trimmed.paired.fastq";
    
    private static final String         OUTFILE_EXTENSION           = ".trimmed.paired.fastq";

    

    private             String          pathToBSMap                 = "";
    private             int             seedSize                    = 16;
    private             String          AlignMode                   = "";
    private             double          noOfMismatches              = 0.08;
    private             int             gapSize                     = 0;
    private             Boolean         threeNucMapping             = false;
    private             int             noOfThreads                 = 4;
    private             String          rootDataFolder              = "";
    private             String          ReferenceGenome             = "";
    private             double          kmerCutoffRatio             = 5e-07;
    private             int             reportRepeatsValue          = 1;
    private             int             qualityThresholdTrimValue   = 0;
    private             int             lowQualityFilterValue       = 5;
    private             int             maxInsertSize               = 1000;
    private             int             minInsertSize               = 28;
    private             int             mapFirstNucleotides         = 160;
    private             int             genomeIndexInterval         = 4;
    private     ArrayList<String>       adapterSequences;
    private             Boolean         trimAdapterSequence         = false;
    private             Boolean         includeRefSeq               = false;
    private             Boolean         skipSAMHeader               = false;
    private             Boolean         reportUnmappedReads         = false;
    private             int             startAtThisRead             = 1;
    private             int             endtAtThisRead              = -1;
    private             String          digestionSite               = "";
    private             int             randomSeed                  = 0;
    private             int             mappingStrand               = 0;
    private             String          transitionMap               = "TC";
    private             int             messageLevel                = 1;
     

    public StepBSMapReads(){
        classSubtype = NGSStepSubclass.DATABLE;
    }
    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepBSMapReads(InputDataForStep sid){
        classSubtype = NGSStepSubclass.DATABLE;
       stepInputData = sid;
    }
    


    @Override
    public String shortStepDescription(){
      return "Perform mapping of BS reads using BSMap";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "Perform mapping of BS reads using BSMap.\n\n"
              + "the calculation of methylation ratios is performed in a separate step.\n\n"
              + "Input is a FASTQ file, raw or trimmed\n"
              + "Output is a alignment file\n";
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
        if(configData.get(ID_DIGESTIONSITE)==null) {
            logger.info("<" + configData.get(ID_DIGESTIONSITE) + "> : Missing Definition in Configuration File");
            logger.error("<" + configData.get(ID_DIGESTIONSITE) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_DIGESTIONSITE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_TRIM_ADAPTERSEQ)==null) {
            logger.info("<" + configData.get(ID_TRIM_ADAPTERSEQ) + "> : Missing Definition in Configuration File");
            logger.error("<" + configData.get(ID_TRIM_ADAPTERSEQ) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_TRIM_ADAPTERSEQ) + "> : Missing Definition in Configuration File");
        }
        
        this.setReferenceGenome((String) configData.get(ID_HOST));
        if(this.getReferenceGenome().length() !=3 ){
            logger.error(ID_HOST + " <" + configData.get(ID_HOST) + "> must be a 3 letter string");            
            throw new IllegalArgumentException(ID_HOST + " <" + configData.get(ID_HOST) + "> must be a 3 letter string");            
        }
        this.setPathToBSMap((String) configData.get(ID_SOFTWARE));
        

        String chk;

        chk = checkParameter("Integer", ID_SEED_SIZE, Integer.toString((Integer)configData.get(ID_SEED_SIZE)), "10", "16", logger);
        if(chk!=null)
            this.setSeedSize((Integer)configData.get(ID_SEED_SIZE));
            
        chk = checkParameter("Double", ID_MISMATCHES, Double.toString((Double)configData.get(ID_MISMATCHES)), "0", "15", logger);
        if(chk!=null)
            this.setNoOfMismatches((Double)configData.get(ID_MISMATCHES));
            
        chk = checkParameter("Integer", ID_GAPSIZE, Integer.toString((Integer)configData.get(ID_GAPSIZE)), "0", "NA", logger);
        if(chk!=null)
            this.setGapSize((Integer)configData.get(ID_GAPSIZE));
                        
        chk = checkParameter("Boolean", ID_THREEMAPPING, Boolean.toString((Boolean)configData.get(ID_THREEMAPPING)), "0", "NA", logger);
        if(chk!=null)
            this.setThreeNucMapping((Boolean)configData.get(ID_THREEMAPPING));
                        
        chk = checkParameter("Integer", ID_THREADS, Integer.toString((Integer)configData.get(ID_THREADS)), "0", "12", logger);
        if(chk!=null)
            this.setNoOfThreads((Integer)configData.get(ID_THREADS));
            
        chk = checkParameter("Double", ID_KMERCUTOFF, Double.toString((Double)configData.get(ID_KMERCUTOFF)), "0", "12", logger);
        if(chk!=null)
            this.setKmerCutoffRatio((Double)configData.get(ID_KMERCUTOFF));
            
        chk = checkParameter("Integer", ID_REPORT_REPEATS, Integer.toString((Integer)configData.get(ID_REPORT_REPEATS)), "0", "2", logger);
        if(chk!=null)
            this.setReportRepeatsValue((Integer)configData.get(ID_REPORT_REPEATS));
            
        chk = checkParameter("Integer", ID_QUAL_THRESHOLD_VALUE, Integer.toString((Integer)configData.get(ID_QUAL_THRESHOLD_VALUE)), "0", "40", logger);
        if(chk!=null)
            this.setQualityThresholdTrimValue((Integer)configData.get(ID_QUAL_THRESHOLD_VALUE));
            
        chk = checkParameter("Integer", ID_LOW_QUAL_VALUE, Integer.toString((Integer)configData.get(ID_LOW_QUAL_VALUE)), "0", "NA", logger);
        if(chk!=null)
            this.setLowQualityFilterValue((Integer)configData.get(ID_LOW_QUAL_VALUE));
            
        chk = checkParameter("Integer", ID_MAX_INSERT_SIZE, Integer.toString((Integer)configData.get(ID_MAX_INSERT_SIZE)), "0", "NA", logger);
        if(chk!=null)
            this.setMaxInsertSize((Integer)configData.get(ID_MAX_INSERT_SIZE));
            
        chk = checkParameter("Integer", ID_MIN_INSERT_SIZE, Integer.toString((Integer)configData.get(ID_MIN_INSERT_SIZE)), "0", "NA", logger);
        if(chk!=null)
            this.setMinInsertSize((Integer)configData.get(ID_MIN_INSERT_SIZE));
            
        chk = checkParameter("Integer", ID_MAP_FIRST_N_NUCS, Integer.toString((Integer)configData.get(ID_MAP_FIRST_N_NUCS)), "0", "NA", logger);
        if(chk!=null)
            this.setMapFirstNucleotides((Integer)configData.get(ID_MAP_FIRST_N_NUCS));
            
        chk = checkParameter("Integer", ID_GENOME_INDEX_INTERVAL, Integer.toString((Integer)configData.get(ID_GENOME_INDEX_INTERVAL)), "1", "16", logger);
        if(chk!=null)
            this.setGenomeIndexInterval((Integer)configData.get(ID_GENOME_INDEX_INTERVAL));

        chk = checkParameter("Boolean", ID_TRIM_ADAPTERSEQ, Boolean.toString((Boolean)configData.get(ID_TRIM_ADAPTERSEQ)), "NA", "NA", logger);
        if(chk!=null){
            this.setTrimAdapterSequence((Boolean)configData.get(ID_TRIM_ADAPTERSEQ));
        }
        if(this.getTrimAdapterSequence()==true){
            this.setAdapterSequences((ArrayList<String>)configData.get(ID_ADAPTERSEQ));
            this.checkAdaptorSequences();
        }

        chk = checkParameter("Boolean", ID_INCLUDE_REF_SEQ, Boolean.toString((Boolean)configData.get(ID_INCLUDE_REF_SEQ)), "NA", "NA", logger);
        if(chk!=null  && Boolean.parseBoolean(chk)==true)
            this.setIncludeRefSeq((Boolean)configData.get(ID_INCLUDE_REF_SEQ));
        
        chk = checkParameter("Boolean", ID_SKIP_SAM_HEADER, Boolean.toString((Boolean)configData.get(ID_SKIP_SAM_HEADER)), "NA", "NA", logger);
        if(chk!=null && Boolean.parseBoolean(chk)==true)
            this.setSkipSAMHeader((Boolean)configData.get(ID_SKIP_SAM_HEADER));
        
        chk = checkParameter("Boolean", ID_REPORT_UNMAPPED_READS, Boolean.toString((Boolean)configData.get(ID_REPORT_UNMAPPED_READS)), "NA", "NA", logger);
        if(chk!=null && Boolean.parseBoolean(chk)==true)
            this.setReportUnmappedReads((Boolean)configData.get(ID_REPORT_UNMAPPED_READS));       
        
        chk = checkParameter("Integer", ID_START_AT_READ_N, Integer.toString((Integer)configData.get(ID_START_AT_READ_N)), "0", "NA", logger);
        if(chk!=null)
            this.setStartAtThisRead((Integer)configData.get(ID_START_AT_READ_N));
            
        chk = checkParameter("Integer", ID_END_AT_READ_N, Integer.toString((Integer)configData.get(ID_END_AT_READ_N)), "-1", "NA", logger);
        if(chk!=null)
            this.setEndtAtThisRead((Integer)configData.get(ID_END_AT_READ_N));
        
        if(this.checkDigestionSite((String)configData.get(ID_DIGESTIONSITE))){
            this.setDigestionSite((String)configData.get(ID_DIGESTIONSITE));
        }

        checkParameter("Integer", ID_RANDOMSEED, Integer.toString((Integer)configData.get(ID_RANDOMSEED)), "0", "NA", logger);
        if(chk!=null)
            this.setRandomSeed((Integer)configData.get(ID_RANDOMSEED));
            
        checkParameter("Integer", ID_MAPPINGSTRAND, Integer.toString((Integer)configData.get(ID_MAPPINGSTRAND)), "0", "NA", logger);
        if(chk!=null)
            this.setMappingStrand((Integer)configData.get(ID_MAPPINGSTRAND));
            
        if(this.checkTransitionMap((String)configData.get(ID_TRANSITIONMAP))){
            this.setTransitionMap((String)configData.get(ID_TRANSITIONMAP));
        }

        checkParameter("Integer", ID_MESSAGELEVEL, Integer.toString((Integer)configData.get(ID_MESSAGELEVEL)), "0", "NA", logger);
        if(chk!=null)
            this.setMessageLevel((Integer)configData.get(ID_MESSAGELEVEL));
            

        logger.info("passed");
    }
    
    
    
    /**
     * check the adaptor sequence is comprised of standard DNA nucleotides
     * 
     * @param seq
     * @return 
     */
    private Boolean checkAdaptorSequences() throws Exception{
        
        for(String adapterSeq: this.getAdapterSequences()){
            for (char ch: adapterSeq.toCharArray()) {
                if("acgtACGT".indexOf(ch) == -1){
                    logger.info("Adaptor Sequence <" + adapterSeq + ">" + " can only contain standard nucleotides (ACGT)");            
                    logger.error("Adaptor Sequence <" + adapterSeq + ">" + " can only contain standard nucleotides (ACGT)");   
                    throw new Exception("Adaptor Sequence <" + adapterSeq + ">" + " can only contain standard nucleotides (ACGT)");

                }
            }
            
        }
        return true;
    }
    



    /**
     * Digestion site should contain standard DNA nucleotides and a '-'
     * to mark the digestion location
     * 
     * @param seq
     * @return 
     */
    private Boolean checkDigestionSite(String seq) throws Exception{
        if(seq.contains("-")==false) return false;
        
        String rawSeq = "";
        for (char ch: seq.toCharArray()) {
            if("acgtACGT-".indexOf(ch) == -1){
                logger.info("Digestion Site <" + seq + ">" + " can only contain standard nucleotides (ACGT)");            
                logger.error("Digestion Site <" + seq + ">" + " can only contain standard nucleotides (ACGT)");            
                throw new Exception("Digestion Site <" + seq + ">" + " can only contain standard nucleotides (ACGT)");
            }
            rawSeq = rawSeq + ch;
        }
        
        if(SimpleSeq.complement(seq.replace("-","")).equals(seq.replace("-",""))==false){
            logger.info("Digestion Site <" + seq + ">" 
                    + " must be palindromic: Complement<" + seq + ">=" + SimpleSeq.complement(seq));            
            logger.error("Digestion Site <" + seq + ">" 
                    + " must be palindromic: Complement<" + seq + ">=" + SimpleSeq.complement(seq));            
            throw new Exception("Digestion Site <" + seq + ">" 
                    + " must be palindromic: Complement<" + seq + ">=" + SimpleSeq.complement(seq));            
        }
        
        return true;
    }
    



    /**
     * Transition Map should comprise a DNA nucleotide pair
     * 
     * @param seq
     * @return 
     */
    private Boolean checkTransitionMap(String seq) throws Exception{
        if(seq.length()!=2) {
            logger.info("Transition Map can only be two nucleotides, <" + seq + "> is "  + seq.length() + " nt");            
            logger.error("Transition Map can only be two nucleotides, <" + seq + "> is "  + seq.length() + " nt");      
            throw new Exception("Transition Map can only be two nucleotides, <" + seq + "> is "  + seq.length() + " nt");
        }
        for (char ch: seq.toCharArray()) {
            if("acgtACGT".indexOf(ch) == -1) {
                logger.info("Transition Map <" + seq + ">" + " can only contain standard nucleotides (ACGT)");            
                logger.error("Transition Map <" + seq + ">" + " can only contain standard nucleotides (ACGT)");            
                throw new Exception("Transition Map <" + seq + ">" + " can only contain standard nucleotides (ACGT)");
            }
        }
        
        return true;
    }
    



    
    /**
     * BSMap single or paired end reads
     * 
     * 
     * @throws IOException 
     */
    @Override
    public void execute() throws IOException{
        /* example call.:
            ./bsmap \
            -a ../../../MetylationDataEPIRA/PilotDNAextractionRnBeads/160125_M01334/1-RA5111B-4nQ_S1_L001_R1_001.paired.fastq.gz \
            -b ../../../MetylationDataEPIRA/PilotDNAextractionRnBeads/160125_M01334/1-RA5111B-4nQ_S1_L001_R2_001.paired.fastq.gz \
            -d ../../../../../../Volumes/LaCie_kari/ref/Homo_sapiens_assembly19.fasta \
            -D C-CGG \
            -p 2 \
            -w 100 \
            -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG \
            -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
            -v 2 \
            -n 1 \
            -o testOut29febSAMformat.sam &
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
                
                String cmdBSMap = "";   
                ArrayList<String> cmd = new ArrayList<>();
                cmd.add(this.getPathToBSMap());
                
                if(this.trimAdapterSequence){
                    fastqFile1out = this.cleanPath(outFolder 
                            + FILESEPARATOR + sampleData.getFastqFile1());
                    fastqFile1in = this.cleanPath(inFolder 
                            + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION_TRIM));                
                }else{
                    fastqFile1out = this.cleanPath(outFolder 
                            + FILESEPARATOR + sampleData.getFastqFile1());
                    fastqFile1in = this.cleanPath(inFolder 
                            + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION_NOTRIM));                
                }
                cmd.add("-a " + fastqFile1in);
                if (sampleData.getFastqFile2() != null){
                    if(this.trimAdapterSequence){
                        fastqFile2out = this.cleanPath(outFolder 
                                + FILESEPARATOR + sampleData.getFastqFile2());
                        fastqFile2in = this.cleanPath(inFolder 
                                + FILESEPARATOR + sampleData.getFastqFile2().replace(".fastq", INFILE_EXTENSION_TRIM));                
                    }else{
                        fastqFile2out = this.cleanPath(outFolder 
                                + FILESEPARATOR + sampleData.getFastqFile2());
                        fastqFile2in = this.cleanPath(inFolder 
                                + FILESEPARATOR + sampleData.getFastqFile2().replace(".fastq", INFILE_EXTENSION_NOTRIM));                
                    }
                    cmd.add("-b " + fastqFile2in);
                }
                
                String pathToGenomeFA = this.cleanPath(getStepInputData().getDataLocations().getGenomeRootFolder()
                        + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_WHOLE_GENSEQ_FA);
                cmd.add("-d " + pathToGenomeFA);
                    
                /**
                 * not sure about the best way to handle this.
                 * It seems that the files that come of the NSC sequences have a general format
                 * SAMPLENAME_SEQUENCERDATA1_SEQUENCEDATA_..., so we split at the first "_"
                 * 
                 */
                String samAlnFile = this.cleanPath(outFolder + FILESEPARATOR 
                        + sampleData.getFastqFile1().split("_")[0].trim()+".sam");
                cmd.add("-o " + samAlnFile);
                
                cmd.add("-s " + this.getSeedSize());
                
                /*  
                  if this value is between 0 and 1, it's interpreted as the mismatch rate w.r.t to the read length.
                  otherwise it's interpreted as the maximum number of mismatches allowed on a read, <=15.
                  example: 
                    -v 5 (max #mismatches = 5), 
                    -v 0.1 (max #mismatches = read_length * 10%)
                    default=0.08.
                */
                cmd.add("-v " + this.getNoOfMismatches());
                
                /*
                  using 3-nucleotide mapping approach. (default: off)
                */
                if(this.getThreeNucMapping())
                    cmd.add("-3");
                
                        
                /*
                  gap size, BSMAP only allow 1 continuous gap (insertion or deletion) 
                  with up to 3 nucleotides. default=0
                  gaps will not be allowed within 6nt of the read edges.
                  the number of mismatches of gapped algnment is calculated as #gap_size+#mismatches+1
                */
                cmd.add("-g " + this.getGapSize());
                
                /*
                  set the cut-off ratio for over-represented kmers, default=1e-06
                  e.g.: -k 1e-6 means the top 0.0001% over-represented kmer will be skipped in alignment
                */
                cmd.add("-k " + this.getKmerCutoffRatio());
                
                /*
                  [0,1,2] how to report repeat hits, 0=none(unique hit/pair only); 
                  1=random one; 2=all(large output file size), default=1.
                */                
                cmd.add("-r " + this.getReportRepeatsValue());

                /*
                  quality threshold in trimming 3'end of reads, 0-40, default=0. (no trim)
                */
                cmd.add("-q " + this.getQualityThresholdTrimValue());
                
                /*
                  filter low-quality reads containing >n Ns, default=5
                */
                cmd.add("-f " + this.getLowQualityFilterValue());
                
                cmd.add("-p " + this.getNoOfThreads());
                cmd.add("-x " + this.getMaxInsertSize());
                cmd.add("-m " + this.getMinInsertSize());
                
                /*
                  map the first N nucleotide of the read, default: 0 (map the whole read).
                */
                cmd.add("-L " + this.getMapFirstNucleotides());

                /*
                  index interval (1~16), meaning the reference genome will be indexed every Nbp, default=4. (WGBS mode)
                  For RRBS mode, index_interval is fixed to 1bp and this command line option is neglected.
                 larger index interval uses memory, and slightly reduces mapping sensitivity. (~0.5% difference)
                 for human genome, -I 16 uses ~5GB, compared with ~9GB at the default -I 4.                
                */
                cmd.add("-I " + this.getGenomeIndexInterval());
                
                if(this.getTrimAdapterSequence()){
                    for(String seq: this.getAdapterSequences()){
                        cmd.add("-A " + seq);
                    }
                }
                
                if(this.getIncludeRefSeq())
                    cmd.add("-R");
                if(!this.getSkipSAMHeader())
                    cmd.add("-H");
                if(this.getReportUnmappedReads())
                    cmd.add("-u");
                
                cmd.add("-B " + this.getStartAtThisRead());
                cmd.add("-E " + this.getEndtAtThisRead());
                cmd.add("-D " + this.getDigestionSite());
                cmd.add("-S " + this.getRandomSeed());
                
                /*
                    -n 0: only map to 2 forward strands, i.e. BSW(++) and BSC(-+) ("Lister protocol")
                    for PE sequencing, map read#1 to ++ and -+, read#2 to +- and --.
                    -n 1: map SE or PE reads to all 4 strands, i.e. ++, +-, -+, -- ("Cokus protocol")
                    default: -n 0. Most bisulfite sequencing data is generated only from forward strands.
                */
                cmd.add("-n " + this.getMappingStrand());
                
                /*
                  set the alignment information for the additional nucleotide transition. <str> is in the form of two different nucleotides,
                  the first one in the reads could be mapped to the second one in the reference sequences.
                  default: -M TC, corresponds to C=>U(T) transition in bisulfite conversion.
                  example: -M GA could be used to detect to A=>I(G) transition in RNA editing.
                */
                cmd.add("-M " + this.getTransitionMap());
                
                /*
                verbose level:  0=no message displayed (quiet mode); 
                                1=major message (default); 
                                2=detailed message                
                */
                cmd.add("-V " + this.getMessageLevel());
                
                /*
                  still not sure how to handle these two parameters
                    -w  <int>   max number of equal best hits to count, smaller will be faster, default=MAXHITS in makefile
                    -z  <int>   base quality, default=33 [Illumina is using 64, Sanger Institute is using 33]
                */                
                                
                cmdBSMap = this.cleanPath(StringUtils.join(cmd, " "));
                logger.info("BSMap command is: " + cmdBSMap);
                Runtime rt = Runtime.getRuntime();
                Process proc = rt.exec(cmdBSMap);
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
                
                
            }
            catch(IOException | InterruptedException ex ){
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
        
        if(new File(this.getPathToBSMap()).exists() == false){
            logger.info("BSMap software not found at location < " + this.getPathToBSMap() +">");
            logger.error("BSMap software not found at location < " + this.getPathToBSMap() +">");
            throw new IOException("BSMap software not found at location < " + this.getPathToBSMap() +">");
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
            
            if(this.trimAdapterSequence){
                //fastqFile1out = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1());
                fastqFile1in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION_TRIM));                
                if (fastqFile1in.toUpperCase().endsWith(INFILE_EXTENSION_TRIM.toUpperCase())==false)
                {
                    logger.info(STEP_ID_STRING + ": incorrect file extension for input file <" 
                      + fastqFile1in + ">. should have <" + INFILE_EXTENSION_TRIM + "> as extension");
                    throw new IOException(STEP_ID_STRING + ": incorrect file extension for input file <" 
                      + fastqFile1in + ">.  \n" 
                      + "should have <" + INFILE_EXTENSION_TRIM + "> as extension");
                }
            }else{
                //fastqFile1out = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1());
                fastqFile1in = this.cleanPath(inFolder 
                        + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION_NOTRIM));                
                if (fastqFile1in.toUpperCase().endsWith(INFILE_EXTENSION_TRIM.toUpperCase())==false)
                {
                    logger.info(STEP_ID_STRING + ": incorrect file extension for input file <" 
                      + fastqFile1in + ">. should have <" + INFILE_EXTENSION_NOTRIM + "> as extension");
                    throw new IOException(STEP_ID_STRING + ": incorrect file extension for input file <" 
                      + fastqFile1in + ">.  \n" 
                      + "should have <" + INFILE_EXTENSION_NOTRIM + "> as extension");
                }
            }
            if (new File(fastqFile1in).exists()==false){
                logger.error(STEP_ID_STRING + ": fastq 1 File <" + fastqFile1in + "> does not exist");
                throw new IOException(STEP_ID_STRING + ": fastq 1 File <" + fastqFile1in + "> does not exist");
            }
            
            
            //Fastq 2
            if (sampleData.getFastqFile2()==null) continue;
            String fastqFile2in = "";
            String fastqFile2out = "";
            if(this.trimAdapterSequence){
                //fastqFile2out = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile2());
                fastqFile2in = this.cleanPath(inFolder 
                        + FILESEPARATOR + sampleData.getFastqFile2().replace(".fastq", INFILE_EXTENSION_TRIM));                
                if (fastqFile2in.toUpperCase().endsWith(INFILE_EXTENSION_TRIM.toUpperCase())==false)
                {
                    logger.error(STEP_ID_STRING + ": incorrect file extension for fastq file 2 <" 
                      + fastqFile2in + ">. should have <" + INFILE_EXTENSION_TRIM + "> as extension");
                    throw new IOException(STEP_ID_STRING + ": incorrect file extension for fastq file 2 <" 
                      + fastqFile2in + ">. \n" 
                      + "should have <" + INFILE_EXTENSION_TRIM + "> as extension");
                }
            }else{
                //fastqFile2out = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile2());
                fastqFile2in = this.cleanPath(inFolder 
                        + FILESEPARATOR + sampleData.getFastqFile2().replace(".fastq", INFILE_EXTENSION_NOTRIM));                
                if (fastqFile2in.toUpperCase().endsWith(INFILE_EXTENSION_NOTRIM.toUpperCase())==false)
                {
                    logger.error(STEP_ID_STRING + ": incorrect file extension for fastq file 2 <" 
                      + fastqFile2in + ">. should have <" + INFILE_EXTENSION_NOTRIM + "> as extension");
                    throw new IOException(STEP_ID_STRING + ": incorrect file extension for fastq file 2 <" 
                      + fastqFile2in + ">. \n" 
                      + "should have <" + INFILE_EXTENSION_NOTRIM + "> as extension");
                }
            }

            
            if ((new File(fastqFile2in)).exists()==false ){
                logger.error(STEP_ID_STRING + ": fastq 2 File <" + fastqFile2in + "> does not exist");
                throw new IOException(STEP_ID_STRING + ": fastq 2 File <" + fastqFile2in + "> does not exist");
            }
                        
            logger.info("passed");
        }

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
        
        //HashMap<String, Object> configData = new HashMap();
        HashMap configData = new HashMap();
        
        configData.put(ID_SOFTWARE, "/usr/local/bin/bsmap");
        configData.put(ID_HOST, "hsa");
        configData.put(ID_SEED_SIZE, 16);
        
        configData.put(ID_MISMATCHES, 0.08);
        configData.put(ID_GAPSIZE, 0);
        configData.put(ID_THREEMAPPING, false);
        configData.put(ID_THREADS, 4);
        configData.put(ID_KMERCUTOFF, 1e-06);
        configData.put(ID_REPORT_REPEATS, 1);
        configData.put(ID_QUAL_THRESHOLD_VALUE, 0);
        configData.put(ID_LOW_QUAL_VALUE, 5);
        configData.put(ID_MAX_INSERT_SIZE, 500);
        configData.put(ID_MIN_INSERT_SIZE, 28);
        configData.put(ID_MAP_FIRST_N_NUCS, 0);
        configData.put(ID_GENOME_INDEX_INTERVAL, 4);
        configData.put(ID_TRIM_ADAPTERSEQ, true);
        configData.put(ID_ADAPTERSEQ, new ArrayList<>(Arrays.asList(
                "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG", 
                "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT")));        
        configData.put(ID_INCLUDE_REF_SEQ, false);
        configData.put(ID_SKIP_SAM_HEADER, false);
        configData.put(ID_REPORT_UNMAPPED_READS, false);
        configData.put(ID_START_AT_READ_N, 1);
        configData.put(ID_END_AT_READ_N, -1);
        configData.put(ID_DIGESTIONSITE, "C-CGG");
        configData.put(ID_RANDOMSEED, 0);
        configData.put(ID_MAPPINGSTRAND, 0);
        configData.put(ID_TRANSITIONMAP, "TC");
        configData.put(ID_MESSAGELEVEL, 1);

        return configData;
        
    }

    
    
    
    
    

    /**
     * @return the pathToBSMap
     */
    public String getPathToBSMap() {
        return pathToBSMap;
    }

    /**
     * @param pathToBSMap the pathToBSMap to set
     */
    public void setPathToBSMap(String pathToBSMap) {
        this.pathToBSMap = pathToBSMap;
    }

    /**
     * @return the seedSize
     */
    public int getSeedSize() {
        return seedSize;
    }

    /**
     * @param seedSize the seedSize to set
     */
    public void setSeedSize(int seedSize) {
        this.seedSize = seedSize;
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
     * @return the noOfMismatches
     */
    public double getNoOfMismatches() {
        return noOfMismatches;
    }

    /**
     * @param noOfMismatches the noOfMismatches to set
     */
    public void setNoOfMismatches(double noOfMismatches) {
        this.noOfMismatches = noOfMismatches;
    }

    /**
     * @return the threeNucMapping
     */
    public Boolean getThreeNucMapping() {
        return threeNucMapping;
    }

    /**
     * @param threeNucMapping the threeNucMapping to set
     */
    public void setThreeNucMapping(Boolean threeNucMapping) {
        this.threeNucMapping = threeNucMapping;
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
     * @return the kmerCutoffRatio
     */
    public double getKmerCutoffRatio() {
        return kmerCutoffRatio;
    }

    /**
     * @param kmerCutoffRatio the kmerCutoffRatio to set
     */
    public void setKmerCutoffRatio(double kmerCutoffRatio) {
        this.kmerCutoffRatio = kmerCutoffRatio;
    }

    /**
     * @return the reportRepeatsValue
     */
    public int getReportRepeatsValue() {
        return reportRepeatsValue;
    }

    /**
     * @param reportRepeatsValue the reportRepeatsValue to set
     */
    public void setReportRepeatsValue(int reportRepeatsValue) {
        this.reportRepeatsValue = reportRepeatsValue;
    }

    /**
     * @return the qualityThresholdTrimValue
     */
    public int getQualityThresholdTrimValue() {
        return qualityThresholdTrimValue;
    }

    /**
     * @param qualityThresholdTrimValue the qualityThresholdTrimValue to set
     */
    public void setQualityThresholdTrimValue(int qualityThresholdTrimValue) {
        this.qualityThresholdTrimValue = qualityThresholdTrimValue;
    }

    /**
     * @return the lowQualityFilterValue
     */
    public int getLowQualityFilterValue() {
        return lowQualityFilterValue;
    }

    /**
     * @param lowQualityFilterValue the lowQualityFilterValue to set
     */
    public void setLowQualityFilterValue(int lowQualityFilterValue) {
        this.lowQualityFilterValue = lowQualityFilterValue;
    }

    /**
     * @return the maxInsertSize
     */
    public int getMaxInsertSize() {
        return maxInsertSize;
    }

    /**
     * @param maxInsertSize the maxInsertSize to set
     */
    public void setMaxInsertSize(int maxInsertSize) {
        this.maxInsertSize = maxInsertSize;
    }

    /**
     * @return the minInsertSize
     */
    public int getMinInsertSize() {
        return minInsertSize;
    }

    /**
     * @param minInsertSize the minInsertSize to set
     */
    public void setMinInsertSize(int minInsertSize) {
        this.minInsertSize = minInsertSize;
    }

    /**
     * @return the mapFirstNucleotides
     */
    public int getMapFirstNucleotides() {
        return mapFirstNucleotides;
    }

    /**
     * @param mapFirstNucleotides the mapFirstNucleotides to set
     */
    public void setMapFirstNucleotides(int mapFirstNucleotides) {
        this.mapFirstNucleotides = mapFirstNucleotides;
    }

    /**
     * @return the genomeIndexInterval
     */
    public int getGenomeIndexInterval() {
        return genomeIndexInterval;
    }

    /**
     * @param genomeIndexInterval the genomeIndexInterval to set
     */
    public void setGenomeIndexInterval(int genomeIndexInterval) {
        this.genomeIndexInterval = genomeIndexInterval;
    }

    /**
     * @return the includeRefSeq
     */
    public Boolean getIncludeRefSeq() {
        return includeRefSeq;
    }

    /**
     * @param includeRefSeq the includeRefSeq to set
     */
    public void setIncludeRefSeq(Boolean includeRefSeq) {
        this.includeRefSeq = includeRefSeq;
    }

    /**
     * @return the skipSAMHeader
     */
    public Boolean getSkipSAMHeader() {
        return skipSAMHeader;
    }

    /**
     * @param skipSAMHeader the skipSAMHeader to set
     */
    public void setSkipSAMHeader(Boolean skipSAMHeader) {
        this.skipSAMHeader = skipSAMHeader;
    }

    /**
     * @return the reportUnmappedReads
     */
    public Boolean getReportUnmappedReads() {
        return reportUnmappedReads;
    }

    /**
     * @param reportUnmappedReads the reportUnmappedReads to set
     */
    public void setReportUnmappedReads(Boolean reportUnmappedReads) {
        this.reportUnmappedReads = reportUnmappedReads;
    }

    /**
     * @return the startAtThisRead
     */
    public int getStartAtThisRead() {
        return startAtThisRead;
    }

    /**
     * @param startAtThisRead the startAtThisRead to set
     */
    public void setStartAtThisRead(int startAtThisRead) {
        this.startAtThisRead = startAtThisRead;
    }

    /**
     * @return the endtAtThisRead
     */
    public int getEndtAtThisRead() {
        return endtAtThisRead;
    }

    /**
     * @param endtAtThisRead the endtAtThisRead to set
     */
    public void setEndtAtThisRead(int endtAtThisRead) {
        this.endtAtThisRead = endtAtThisRead;
    }

    /**
     * @return the digestionSite
     */
    public String getDigestionSite() {
        return digestionSite;
    }

    /**
     * @param digestionSite the digestionSite to set
     */
    public void setDigestionSite(String digestionSite) {
        this.digestionSite = digestionSite;
    }

    /**
     * @return the randomSeed
     */
    public int getRandomSeed() {
        return randomSeed;
    }

    /**
     * @param randomSeed the randomSeed to set
     */
    public void setRandomSeed(int randomSeed) {
        this.randomSeed = randomSeed;
    }

    /**
     * @return the mappingStrand
     */
    public int getMappingStrand() {
        return mappingStrand;
    }

    /**
     * @param mappingStrand the mappingStrand to set
     */
    public void setMappingStrand(int mappingStrand) {
        this.mappingStrand = mappingStrand;
    }

    /**
     * @return the transitionMap
     */
    public String getTransitionMap() {
        return transitionMap;
    }

    /**
     * @param transitionMap the transitionMap to set
     */
    public void setTransitionMap(String transitionMap) {
        this.transitionMap = transitionMap;
    }

    /**
     * @return the gapSize
     */
    public int getGapSize() {
        return gapSize;
    }

    /**
     * @param gapSize the gapSize to set
     */
    public void setGapSize(int gapSize) {
        this.gapSize = gapSize;
    }

    /**
     * @return the messageLevel
     */
    public int getMessageLevel() {
        return messageLevel;
    }

    /**
     * @param messageLevel the messageLevel to set
     */
    public void setMessageLevel(int messageLevel) {
        this.messageLevel = messageLevel;
    }

    /**
     * @return the trimAdapterSequence
     */
    public Boolean getTrimAdapterSequence() {
        return trimAdapterSequence;
    }

    /**
     * @param trimAdapterSequence the trimAdapterSequence to set
     */
    public void setTrimAdapterSequence(Boolean trimAdapterSequence) {
        this.trimAdapterSequence = trimAdapterSequence;
    }

    /**
     * @return the adapterSequences
     */
    public ArrayList<String> getAdapterSequences() {
        return adapterSequences;
    }

    /**
     * @param adapterSequences the adapterSequences to set
     */
    public void setAdapterSequences(ArrayList<String> adapterSequences) {
        this.adapterSequences = adapterSequences;
    }


    
    
    
    
    
}
