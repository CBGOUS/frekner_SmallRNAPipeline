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
 * HISAT2 mapping for RNA Seq
 *
 * The step will only accept FASTQ format for sequence files
 * all parameters specified in the user manual may be applied here
 * for parameters that are flags, these may be turned off by setting 'NO' in the
 * YAML file.
 * 
 * The step is only written to handle Illumina reads
 * This step handles both single and paired end reads
 *
 * @author sr
 */
public class StepHisat2MapReads extends NGSStep{

    static Logger logger = LogManager.getLogger();

    public static final NGSStepSubclass STEP_SUBCLASS                       = NGSStepSubclass.DATABLE;

    public  static final String     STEP_ID_STRING                          = "HISAT2Mapping";

    // step specific parameters
    private static final String     ID_SOFTWARE                             = "path_to_hisat2";
    private static final String     ID_REF_GENOME                           = "host";
    
    // HISAT2 specific parameters
    
    // Input options
    private static final String     HS2_SKIP_FIRST_NREADS                   = "skip";
    private static final String     HS2_ALIGN_FIRST_NREADS                  = "qupto";
    private static final String     HS2_TRIM_BASES_FROM_5PRIME_END          = "trim5";
    private static final String     HS2_TRIM_BASES_FROM_3PRIME_END          = "trim3";
    
    // Alignment options
    private static final String     HS2_N_CEILING_F                         = "n-ceil-fntype";
    private static final String     HS2_N_CEILING_C                         = "n-ceil-constant";
    private static final String     HS2_N_CEILING_E                         = "n-ceil-coeff";
    private static final String     HS2_NOFW                                = "nofw";
    private static final String     HS2_NOFC                                = "nofc";
    private static final String     HS2_IGNORE_QUALITY_VALUES               = "ignore-quals";
    
    // Scoring options
    private static final String     HS2_PENALTY_MAX_MISMATCH                = "mp-max";
    private static final String     HS2_PENALTY_MIN_MISMATCH                = "mp-min";
    private static final String     HS2_PENALTY_MAX_SOFTCLIP_MISMATCH       = "sp-max";
    private static final String     HS2_PENALTY_MIN_SOFTCLIP_MISMATCH       = "sp-min";
    private static final String     HS2_NO_SOFTCLIP                         = "no-softclip";
    private static final String     HS2_PENALTY_AMBIGUOUS_CHARACTER         = "np";
    private static final String     HS2_PENALTY_READGAP_OPEN                = "readdgap-open";
    private static final String     HS2_PENALTY_READGAP_EXTEND              = "readdgap-extend";
    private static final String     HS2_PENALTY_REFGAP_OPEN                 = "refgap-open";
    private static final String     HS2_PENALTY_REFGAP_EXTEND               = "refgap-extend";
    private static final String     HS2_MIN_ALN_SCORE_FN_C                  = "score-min-fntype";
    private static final String     HS2_MIN_ALN_SCORE_FN_L                  = "score-min-constant";
    private static final String     HS2_MIN_ALN_SCORE_FN_S                  = "score-min-coeff";
    
    // Spliced alignment options
    private static final String     HS2_PENALTY_CANONICAL_SPLICE            = "pen-cansplice";
    private static final String     HS2_PENALTY_NONCANON_SPLICE             = "pen-noncansplice";
    private static final String     HS2_PENALTY_LONGINTRON_CAN_SPLICE_F     = "pen-canintronlen-fntype";
    private static final String     HS2_PENALTY_LONGINTRON_CAN_SPLICE_C     = "pen-canintronlen-constant";
    private static final String     HS2_PENALTY_LONGINTRON_CAN_SPLICE_E     = "pen-canintronlen-coeff";
    private static final String     HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_F  = "pen-noncanintronlen-fntype";
    private static final String     HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_C  = "pen-noncanintronlen-constant";
    private static final String     HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_E  = "pen-noncanintronlen-coeff";
    private static final String     HS2_MIN_INTRON_LEN                      = "min-intronlen";
    private static final String     HS2_MAX_INTRON_LEN                      = "max-intronlen";
    private static final String     HS2_KNOWN_SPLICES_SITES_INFILE          = "known-splicesite-infile";
    private static final String     HS2_NOVEL_SPLICES_SITES_OUTFILE         = "novel-splicesite-outfile";
    private static final String     HS2_NOVEL_SPLICES_SITES_INFILE          = "novel-splicesite-infile";
    private static final String     HS2_NO_TEMP_SPLICESITE                  = "no-temp-splicesite";
    private static final String     HS2_NO_SPLICED_ALIGNMENT                = "no-spliced-alignment";
    private static final String     HS2_RNA_STRANDEDNESS                    = "rna-strandness";
    private static final String     HS2_TRANSCRIPTOME_MAPPING_ONLY          = "transcript-mapping-only";
    private static final String     HS2_DOWNSTREAM_TRANSCRIPTOME_ASSEMBLY   = "downstream-transcriptome-assembly";
    private static final String     HS2_DTA_CUFFLINKS                       = "dta-cufflinks";
    
    // Reporting Options
    private static final String     HS2_NUM_DISTINCT_PRIMARY_ALNS           = "k";
    private static final String     HS2_MAX_NUM_SEEDS                       = "max-seeds";
    private static final String     HS2_REPORT_SECONDARY_ALNS               = "secondary";
    
    // Paired End options
    private static final String     HS2_MIN_FRAG_LEN_FOR_PAIREDEND          = "minins";
    private static final String     HS2_MAX_FRAG_LEN_FOR_PAIREDEND          = "maxins";
    private static final String     HS2_MATE_ORIENTATION_FR                 = "fr";
    private static final String     HS2_MATE_ORIENTATION_RF                 = "rf";
    private static final String     HS2_MATE_ORIENTATION_FF                 = "ff";
    private static final String     HS2_NO_MIXED_ALN                        = "no-mixed";
    private static final String     HS2_NO_DISCORDANT_ALN                   = "no-discordant";
    
    // Output options
    private static final String     HS2_PRINT_WALL_CLOCK                    = "time";
    private static final String     HS2_WRITE_UNPAIRED_READS_NOALN          = "un";
    private static final String     HS2_WRITE_UNPAIRED_READS_NOALN_GZ       = "un-gz";
    private static final String     HS2_WRITE_UNPAIRED_READS_NOALN_BZ2      = "un-bz2";
    private static final String     HS2_WRITE_UNPAIRED_READS_ALN            = "al";
    private static final String     HS2_WRITE_UNPAIRED_READS_ALN_GZ         = "al-gz";
    private static final String     HS2_WRITE_UNPAIRED_READS_ALN_BZ2        = "al-bz2";
    private static final String     HS2_WRITE_UNPAIRED_READS_DISC_ALN       = "un-conc";
    private static final String     HS2_WRITE_UNPAIRED_READS_DISC_ALN_GZ    = "un-conc-gz";
    private static final String     HS2_WRITE_UNPAIRED_READS_DISC_ALN_BZ2   = "un-conc-bz2";
    private static final String     HS2_WRITE_UNPAIRED_READS_CONC_ALN       = "al-conc";
    private static final String     HS2_WRITE_UNPAIRED_READS_CONC_ALN_GZ    = "al-conc-gz";
    private static final String     HS2_WRITE_UNPAIRED_READS_CONC_ALN_BZ2   = "al-conc-bz2";
    
    private static final String     HS2_QUIET_MODE                          = "quiet";
    private static final String     HS2_ALN_SUMMARY_FILE                    = "summary-file";
    private static final String     HS2_ALN_NEW_SUMMARY                     = "new-summary";
    private static final String     HS2_HISAT2_METRICS_FILE                 = "met-file";
    private static final String     HS2_HISAT2_METRICS_STDERR               = "met-stderr";
    private static final String     HS2_HISAT2_METRICS_TIME_INTERVAL        = "met";
    
    // SAM options
    private static final String     HS2_SAM_SUPPRESS_NOALN                  = "no-unal";
    private static final String     HS2_SAM_SUPPRESS_HEADER_LINES           = "no-hd";
    private static final String     HS2_SAM_SUPPRESS_SQ_LINES               = "no-sq";
    private static final String     HS2_SAM_SET_READ_GROUP_ID               = "rg-id";
    private static final String     HS2_SAM_ADD_READ_GROUP_TEXT             = "rg";
    private static final String     HS2_SAM_REM_CHR_TEXT                    = "remove-chrname";
    private static final String     HS2_SAM_ADD_CHR_TEXT                    = "add-chrname";
    private static final String     HS2_SAM_OMIT_SEC_AND_QUAL_STRINGS       = "omit-sec-seq";
    
    // Performance options
    private static final String     HS2_OVERRIDE_OFFRATE                    = "offrate";
    private static final String     HS2_NTHREADS                            = "threads";
    private static final String     HS2_REORDER                             = "reorder";
    private static final String     HS2_MEMORY_MAPPED_IO                    = "mm";   
    
    // Other options
    private static final String     HS2_QC_FILTER                           = "qc_filter";
    private static final String     HS2_SEED                                = "qseed";
    private static final String     HS2_NON_DETERMINISTIC                   = "non-deterministic";
    private static final String     HS2_VERSION                             = "version";
    private static final String     HS2_HELP                                = "help";
    
    // File extensions
    private static final String     INFILE_EXTENSION                        = ".fastq";
    private static final String     UNALN_EXTENSION                         = ".unaln.fastq";
    private static final String     ALN_EXTENSION                           = ".aln.fastq";
    private static final String     DISCONCORD_ALN_EXTENSION                = ".dconcordaln.fastq";
    private static final String     CONCORD_ALN_EXTENSION                   = ".concordaln.fastq";
    private static final String     SAM_GENALN_EXTENSION                    = ".sam";

    // Step specific variables
    private             String      mappingSoftware                 = "";
    private             String      rootDataFolder                  = "";
    private             String      ReferenceGenome                 = "";
    
    // HISAT2 specific variables
    // Input variables
    private             int         skipFirstNReads                 = 0;
    private             int         alignFirstNReads                = -1;
    private             int         trimNBasesFrom5Prime            = 0;
    private             int         trimNBasesFrom3Prime            = 0;
    
    // Alignment variables
    private             String      nCeilingFn                      = "-L";
    private             double      nCeilingConst                   = 0.0;
    private             double      nCeilingCoeff                   = 0.15;
    private             Boolean     nofw                            = false;
    private             Boolean     nofc                            = false;
    private             Boolean     ignoreQuals                     = false;
    
    //Scoring variables
    private             int         maxMismatchPenalty              = 6;
    private             int         minMismatchPenalty              = 2;
    private             int         maxsoftClipPenalty              = 2;
    private             int         minsoftClipPenalty              = 1;
    private             Boolean     allowSoftClip                   = true;
    private             int         penaltyAmbiguousChar            = 1;
    private             int         readGapOpenPenalty              = 5;
    private             int         readGapExtendPenalty            = 3;    
    private             int         refGapOpenPenalty               = 5;
    private             int         refGapExtendPenalty             = 3;
    private             String      scoreMinFun                     = "L";
    private             double      scoreMinConst                   = 0;
    private             double      scoreMinCoeff                   = 0.2;
    
    //Spliced Alignment variables
    private             int         penaltyCanSplice                = 0;
    private             int         penaltyNonCanSplice             = 12;
    private             String      penaltyLongCanSpiceFun          = "G";
    private             double      penaltyLongCanSpiceConst        = -8;
    private             double      penaltyLongCanSpiceCoeff        = 1;
    private             String      penaltyLongNonCanSpiceFun       = "G";
    private             double      penaltyLongNonCanSpiceConst     = -8;
    private             double      penaltyLongNonCanSpiceCoeff     = 1;
    private             int         minIntronLength                 = 20;
    private             int         maxIntronLength                 = 500000;
    private             String      knownSpliceSitesFileIn          = "";
    private             String      novelSpliceSitesFileIn          = "";
    private             String      novelSpliceSitesFileOut         = "";
    private             Boolean     disableSplicedAligned           = false;
    private             String      rnaStrandness                   = "unstranded";
    private             Boolean     transcriptomeAnalysisOnly       = false;
    private             Boolean     dta                             = false;
    private             Boolean     dtaCufflinks                    = false;
    private             Boolean     disableTemplateAdjustment       = false;
    
    //Reporting variables
    private             int         numOfDistinctPrimaryAligns      = 5;
    private             int         maxNumSeeds                     = 5;
    private             Boolean     reportSecondaryAligns           = false;
    
    //Paired End variables
    private             int         minFraglenPairedAlign           = 0;
    private             int         maxFraglenPairedAlign           = 500;
    private             Boolean     mateOrientationFR               = true;
    private             Boolean     mateOrientationRF               = false;
    private             Boolean     mateOrientationFF              = false;
    private             Boolean     noMixedAligns                   = false;
    private             Boolean     discordantAligns                = true;
    
    //Output variables
    private             Boolean     wallClock                       = false;
    private             Boolean     writeUnpairedReads              = true;
    private             Boolean     writeAlignedReads               = true;
    private             Boolean     writeConcordantReads            = true;
    private             Boolean     writeNonconcordantReads         = true;
    private             Boolean     gzipReadFiles                   = true;
    private             Boolean     bz2ReadFiles                    = false;
    
    private             Boolean     quietMode                       = false;
    private             Boolean     writeAlignmentSummary           = true;
    private             Boolean     writeAlignmentNewStyle          = true;
    private             Boolean     writeHisatMetrics               = true;
    private             Boolean     writeHisatMetricsStderr         = false;
    private             int         metWriteInterval                = 1;
    
    //SAM variables
    private             Boolean     samNoUnaligedReads              = false;
    private             Boolean     samNoHeaderLines                = false;
    private             Boolean     samNoSeqHeaderLines             = false;
    private             String      samReadGroupID                  = "";
    private             String      samReadGroupText                = "";
    private             Boolean     samRemoveChrStringFromRef       = false;
    private             Boolean     samAddChrStringFromRef          = false;
    private             Boolean     samOmitSeqAndQualStrin          = false;
    
    //Performance variables
    private             int         indexOffrate                    = -1;
    private             int         numOfThreads                    = 1;
    private             Boolean     reorderSAMRecords               = false;
    private             Boolean     useMemoryMappedIO               = false;
    
    //Other variables
    private             Boolean     qcFiltering                     = false;
    private             int         seed                            = 0;
    private             Boolean     nonDeterministicAlignment       = false;
    
    
    
    
    
    private             String  samAligmentFile                 = "";
    private             String  fastqInputFile                  = "";
    private             String  fastqGenomeAln                  = "";
    private             String  fastqGenomeUnAln                = "";
    
    private  ArrayList<String>  mapGenStdErr;

    public StepHisat2MapReads() {
        classSubtype = NGSStepSubclass.DATABLE;
    }
    
    /**
     *
     * @param sid StepInputData
     *
     */
    public StepHisat2MapReads(InputDataForStep sid) {
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    
    
    @Override
    public String shortStepDescription(){
      return "Performs HISAT2 mapping of paired reads to reference genome";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "Performs HISAT2 mapping of paired reads to refence genome.\n"
              + "The step requires the Bowtie2 software to be installed.\n"
              + "The location can be specified in the YAML configuration file\n"
              + "along with mapping parameters and reference genome. \n";
    }

    
    
    
    
    /**
     * in this method we are simply checking that the configuration file 
     * has all the entries we need. We don't check if the values are acceptable
     * - that is performed in verifyInputData().
     * 
     * here, we are only check that the minimum data required to perform 
     * HISAT2 mapping is specified in the YAML
     * 
     * this comprises:
     * 
     *  - reference genome
     *  - path to software
     * 
     * We also check that any other specified HISAT2 parameters have acceptable values
     * As there are so many parameters, we break them down into separate methods
     * for each of the categories defined in the HISAT2 manual.
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{
        int level = 2;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Input parameters");

        logger.info(indent + STEP_ID_STRING + ": verify configuration data");
        
        if(configData.get(ID_SOFTWARE)==null) {
            throw new NullPointerException("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_REF_GENOME)==null) {
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
                

        this.setReferenceGenome((String) configData.get(ID_REF_GENOME));
        if(this.getReferenceGenome().length() !=3 ){
            throw new IllegalArgumentException(ID_REF_GENOME + " <" + (String) configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
        }
        this.setMappingSoftware((String) configData.get(ID_SOFTWARE));
        

        this.parseHisat2InputConfigurationData(configData);
        this.parseHisat2AlignmentConfigurationData(configData);
        this.parseHisat2ScoringConfigurationData(configData);
        this.parseHisat2SplicedConfigurationData(configData);
        this.parseHisat2ReportingConfigurationData(configData);
        this.parseHisat2PairedEndConfigurationData(configData);
        this.parseHisat2OutputConfigurationData(configData);
        this.parseHisat2SAMConfigurationData(configData);
        this.parseHisat2PerformanceConfigurationData(configData);
        this.parseHisat2OtherConfigurationData(configData);
        
        logger.info(indent + "passed");
    }
    
    /**
     * parse HISAT2 parameters defining input settings
     * 
     * @param configData
     * @throws Exception 
     */
    public void parseHisat2InputConfigurationData(HashMap configData) throws Exception{
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Input parameters");

        String chkInput;
        
        chkInput = checkParameter("Integer", HS2_SKIP_FIRST_NREADS, 
                Integer.toString((Integer)configData.get(HS2_SKIP_FIRST_NREADS)), 
                "1", "NA", logger);
        if(chkInput!=null)
            this.setSkipFirstNReads((Integer)configData.get(HS2_SKIP_FIRST_NREADS));

        chkInput = checkParameter("Integer", HS2_ALIGN_FIRST_NREADS, 
                Integer.toString((Integer)configData.get(HS2_ALIGN_FIRST_NREADS)), 
                "1", "NA", logger);
        if(chkInput!=null)
            this.setSkipFirstNReads((Integer)configData.get(HS2_ALIGN_FIRST_NREADS));

        chkInput = checkParameter("Integer", HS2_TRIM_BASES_FROM_5PRIME_END, 
                Integer.toString((Integer)configData.get(HS2_TRIM_BASES_FROM_5PRIME_END)), 
                "1", "NA", logger);
        if(chkInput!=null)
            this.setTrimNBasesFrom5Prime((Integer)configData.get(HS2_TRIM_BASES_FROM_5PRIME_END));

        chkInput = checkParameter("Integer", HS2_TRIM_BASES_FROM_3PRIME_END, 
                Integer.toString((Integer)configData.get(HS2_TRIM_BASES_FROM_3PRIME_END)), 
                "1", "NA", logger);
        if(chkInput!=null)
            this.setTrimNBasesFrom3Prime((Integer)configData.get(HS2_TRIM_BASES_FROM_3PRIME_END));


        logger.info(indent + "passed");
        
    }
    
    
    /**
     * parse HISAT2 parameters defining input settings
     * 
     * @param configData
     * @throws Exception 
     */
    public void parseHisat2AlignmentConfigurationData(HashMap configData) throws Exception{
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Alignment parameters");
        
        String chkAln;
        
        chkAln = checkParameter("String", HS2_N_CEILING_F, 
                configData.get(HS2_N_CEILING_F).toString(), 
                "1", "1", logger);
        if(chkAln!=null)
            this.setnCeilingFn((String)configData.get(HS2_N_CEILING_F));

        chkAln = checkParameter("Double", HS2_N_CEILING_C, 
                Double.toString((Double)configData.get(HS2_N_CEILING_C)), 
                "1", "NA", logger);
        if(chkAln!=null)
            this.setnCeilingConst((Double)configData.get(HS2_N_CEILING_C));
        
        chkAln = checkParameter("Double", HS2_N_CEILING_E, 
                Double.toString((Double)configData.get(HS2_N_CEILING_E)), 
                "1", "NA", logger);
        if(chkAln!=null)
            this.setnCeilingCoeff((Double)configData.get(HS2_N_CEILING_E));
        
        if(configData.get(HS2_NOFW)!=null)
            this.setNofw(true);
        if(configData.get(HS2_NOFC)!=null)
            this.setNofc(true);
        if(configData.get(HS2_IGNORE_QUALITY_VALUES)!=null)
            this.setIgnoreQuals(true);
                        
        
        logger.info(indent + "passed");
        
    }
    
    
    /**
     * parse HISAT2 parameters defining input settings
     * 
     * @param configData
     * @throws Exception 
     */
    public void parseHisat2ScoringConfigurationData(HashMap configData) throws Exception{
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Scoring parameters");
        
        String chkScoring = "";
        chkScoring = checkParameter("Integer", HS2_PENALTY_MAX_MISMATCH, 
                Integer.toString((Integer)configData.get(HS2_PENALTY_MAX_MISMATCH)), 
                "1", "NA", logger);
        if(chkScoring!=null)
            this.setMaxMismatchPenalty((Integer)configData.get(HS2_PENALTY_MAX_MISMATCH));
        
        chkScoring = checkParameter("Integer", HS2_PENALTY_MIN_MISMATCH, 
                Integer.toString((Integer)configData.get(HS2_PENALTY_MIN_MISMATCH)), 
                "1", "NA", logger);
        if(chkScoring!=null)
            this.setMinMismatchPenalty((Integer)configData.get(HS2_PENALTY_MIN_MISMATCH));
        
        chkScoring = checkParameter("Integer", HS2_PENALTY_MAX_SOFTCLIP_MISMATCH, 
                Integer.toString((Integer)configData.get(HS2_PENALTY_MAX_SOFTCLIP_MISMATCH)), 
                "1", "NA", logger);
        if(chkScoring!=null)
            this.setMaxsoftClipPenalty((Integer)configData.get(HS2_PENALTY_MAX_SOFTCLIP_MISMATCH));
                
        chkScoring = checkParameter("Integer", HS2_PENALTY_MIN_SOFTCLIP_MISMATCH, 
                Integer.toString((Integer)configData.get(HS2_PENALTY_MIN_SOFTCLIP_MISMATCH)), 
                "1", "NA", logger);
        if(chkScoring!=null)
            this.setMinsoftClipPenalty((Integer)configData.get(HS2_PENALTY_MIN_SOFTCLIP_MISMATCH));
        
        if(configData.get(HS2_NO_SOFTCLIP)!=null)
            this.setAllowSoftClip(true);
                
        chkScoring = checkParameter("Integer", HS2_PENALTY_AMBIGUOUS_CHARACTER, 
                Integer.toString((Integer)configData.get(HS2_PENALTY_AMBIGUOUS_CHARACTER)), 
                "1", "NA", logger);
        if(chkScoring!=null)
            this.setPenaltyAmbiguousChar((Integer)configData.get(HS2_PENALTY_AMBIGUOUS_CHARACTER));
        
        chkScoring = checkParameter("Integer", HS2_PENALTY_READGAP_OPEN, 
                Integer.toString((Integer)configData.get(HS2_PENALTY_READGAP_OPEN)), 
                "1", "NA", logger);
        if(chkScoring!=null)
            this.setReadGapOpenPenalty((Integer)configData.get(HS2_PENALTY_READGAP_OPEN));
        
        chkScoring = checkParameter("Integer", HS2_PENALTY_READGAP_EXTEND, 
                Integer.toString((Integer)configData.get(HS2_PENALTY_READGAP_EXTEND)), 
                "1", "NA", logger);
        if(chkScoring!=null)
            this.setReadGapExtendPenalty((Integer)configData.get(HS2_PENALTY_READGAP_EXTEND));
        
        chkScoring = checkParameter("Integer", HS2_PENALTY_REFGAP_OPEN, 
                Integer.toString((Integer)configData.get(HS2_PENALTY_REFGAP_OPEN)), 
                "1", "NA", logger);
        if(chkScoring!=null)
            this.setRefGapOpenPenalty((Integer)configData.get(HS2_PENALTY_REFGAP_OPEN));
        
        chkScoring = checkParameter("Integer", HS2_PENALTY_REFGAP_EXTEND, 
                Integer.toString((Integer)configData.get(HS2_PENALTY_REFGAP_EXTEND)), 
                "1", "NA", logger);
        if(chkScoring!=null)
            this.setRefGapExtendPenalty((Integer)configData.get(HS2_PENALTY_REFGAP_EXTEND));


        
        chkScoring = checkParameter("String", HS2_MIN_ALN_SCORE_FN_C, 
                configData.get(HS2_MIN_ALN_SCORE_FN_C).toString(), 
                "1", "1", logger);
        if(chkScoring!=null)
            this.setScoreMinFun((String)configData.get(HS2_MIN_ALN_SCORE_FN_C));

        chkScoring = checkParameter("Double", HS2_MIN_ALN_SCORE_FN_L, 
                Double.toString((Double)configData.get(HS2_MIN_ALN_SCORE_FN_L)), 
                "1", "NA", logger);
        if(chkScoring!=null)
            this.setScoreMinConst((Double)configData.get(HS2_MIN_ALN_SCORE_FN_L));
        
        chkScoring = checkParameter("Double", HS2_MIN_ALN_SCORE_FN_S, 
                Double.toString((Double)configData.get(HS2_MIN_ALN_SCORE_FN_S)), 
                "1", "NA", logger);
        if(chkScoring!=null)
            this.setScoreMinCoeff((Double)configData.get(HS2_MIN_ALN_SCORE_FN_S));
        
        logger.info(indent + "passed");
    }
    
    
    /**
     * parse HISAT2 parameters defining input settings
     * 
     * @param configData
     * @throws Exception 
     */
    public void parseHisat2SplicedConfigurationData(HashMap configData) throws Exception{
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Spliced parameters");
        
        String chkSplice = "";
        chkSplice = checkParameter("Integer", HS2_PENALTY_CANONICAL_SPLICE, 
                Integer.toString((Integer)configData.get(HS2_PENALTY_CANONICAL_SPLICE)), 
                "1", "NA", logger);
        if(chkSplice!=null)
            this.setPenaltyCanSplice((Integer)configData.get(HS2_PENALTY_CANONICAL_SPLICE));
        
        chkSplice = checkParameter("Integer", HS2_PENALTY_NONCANON_SPLICE, 
                Integer.toString((Integer)configData.get(HS2_PENALTY_NONCANON_SPLICE)), 
                "1", "NA", logger);
        if(chkSplice!=null)
            this.setPenaltyNonCanSplice((Integer)configData.get(HS2_PENALTY_NONCANON_SPLICE));
        

        chkSplice = checkParameter("String", HS2_PENALTY_LONGINTRON_CAN_SPLICE_F, 
                configData.get(HS2_PENALTY_LONGINTRON_CAN_SPLICE_F).toString(), 
                "1", "1", logger);
        if(chkSplice!=null)
            this.setPenaltyLongCanSpiceFun((String)configData.get(HS2_PENALTY_LONGINTRON_CAN_SPLICE_F));

        chkSplice = checkParameter("Double", HS2_PENALTY_LONGINTRON_CAN_SPLICE_C, 
                Double.toString((Double)configData.get(HS2_PENALTY_LONGINTRON_CAN_SPLICE_C)), 
                "1", "NA", logger);
        if(chkSplice!=null)
            this.setPenaltyLongCanSpiceConst((Double)configData.get(HS2_PENALTY_LONGINTRON_CAN_SPLICE_C));
        
        chkSplice = checkParameter("Double", HS2_PENALTY_LONGINTRON_CAN_SPLICE_E, 
                Double.toString((Double)configData.get(HS2_PENALTY_LONGINTRON_CAN_SPLICE_E)), 
                "1", "NA", logger);
        if(chkSplice!=null)
            this.setPenaltyLongCanSpiceCoeff((Double)configData.get(HS2_PENALTY_LONGINTRON_CAN_SPLICE_E));
        
        
        chkSplice = checkParameter("String", HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_F, 
                configData.get(HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_F).toString(), 
                "1", "1", logger);
        if(chkSplice!=null)
            this.setPenaltyLongNonCanSpiceFun((String)configData.get(HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_F));

        chkSplice = checkParameter("Double", HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_C, 
                Double.toString((Double)configData.get(HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_C)), 
                "1", "NA", logger);
        if(chkSplice!=null)
            this.setPenaltyLongNonCanSpiceConst((Double)configData.get(HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_C));
        
        chkSplice = checkParameter("Double", HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_E, 
                Double.toString((Double)configData.get(HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_E)), 
                "1", "NA", logger);
        if(chkSplice!=null)
            this.setPenaltyLongNonCanSpiceCoeff((Double)configData.get(HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_E));
        

        chkSplice = checkParameter("Integer", HS2_MIN_INTRON_LEN, 
                Integer.toString((Integer)configData.get(HS2_MIN_INTRON_LEN)), 
                "1", "NA", logger);
        if(chkSplice!=null)
            this.setMinIntronLength((Integer)configData.get(HS2_MIN_INTRON_LEN));
        
        chkSplice = checkParameter("Integer", HS2_MAX_INTRON_LEN, 
                Integer.toString((Integer)configData.get(HS2_MAX_INTRON_LEN)), 
                "1", "NA", logger);
        if(chkSplice!=null)
            this.setMaxIntronLength((Integer)configData.get(HS2_MAX_INTRON_LEN));

        
        if(configData.get(HS2_NO_SPLICED_ALIGNMENT)!=null)
            this.setDisableSplicedAligned(true);                
        if(configData.get(HS2_TRANSCRIPTOME_MAPPING_ONLY)!=null) 
            this.setTranscriptomeAnalysisOnly(true);                
        if(configData.get(HS2_DOWNSTREAM_TRANSCRIPTOME_ASSEMBLY)!=null) 
            this.setDta(true);
        if(configData.get(HS2_DTA_CUFFLINKS)!=null) 
            this.setDtaCufflinks(true);                
 
        if(configData.get(HS2_KNOWN_SPLICES_SITES_INFILE)!=null)
            this.setKnownSpliceSitesFileIn(((String)configData.get(HS2_KNOWN_SPLICES_SITES_INFILE)).trim());
        if(configData.get(HS2_NOVEL_SPLICES_SITES_INFILE)!=null)
            this.setNovelSpliceSitesFileIn(((String)configData.get(HS2_NOVEL_SPLICES_SITES_INFILE)).trim());
        if(configData.get(HS2_NOVEL_SPLICES_SITES_OUTFILE)!=null)
            this.setNovelSpliceSitesFileOut(((String)configData.get(HS2_NOVEL_SPLICES_SITES_OUTFILE)).trim());

        // the veracity of this value needs to be checked at run time as the allowed 
        // values depend on whether there are single or paired end reads
        if(configData.get(HS2_RNA_STRANDEDNESS)!=null)
            this.setRnaStrandness(((String)configData.get(HS2_RNA_STRANDEDNESS)).trim());
        logger.info(indent + "passed");
        
    }
    
    
    /**
     * parse HISAT2 parameters defining input settings
     * 
     * @param configData
     * @throws Exception 
     */
    public void parseHisat2ReportingConfigurationData(HashMap configData) throws Exception{
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Reporting parameters");
        
        String chkReport = "";
        chkReport = checkParameter("Integer", HS2_NUM_DISTINCT_PRIMARY_ALNS, 
                Integer.toString((Integer)configData.get(HS2_NUM_DISTINCT_PRIMARY_ALNS)), 
                "1", "NA", logger);
        if(chkReport!=null)
            this.setNumOfDistinctPrimaryAligns((Integer)configData.get(HS2_NUM_DISTINCT_PRIMARY_ALNS));

        chkReport = checkParameter("Integer", HS2_MAX_NUM_SEEDS, 
                Integer.toString((Integer)configData.get(HS2_MAX_NUM_SEEDS)), 
                "1", "NA", logger);
        if(chkReport!=null)
            this.setMaxNumSeeds((Integer)configData.get(HS2_MAX_NUM_SEEDS));

        if((Boolean)configData.get(HS2_REPORT_SECONDARY_ALNS)!=null)
            this.setReportSecondaryAligns(true);

        logger.info(indent + "passed");
        
    }
    
    
    /**
     * parse HISAT2 parameters defining input settings
     * 
     * @param configData
     * @throws Exception 
     */
    public void parseHisat2PairedEndConfigurationData(HashMap configData) throws Exception{
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 PairedEnd parameters");
        
        String chkPair = "";
        chkPair = checkParameter("Integer", HS2_MIN_FRAG_LEN_FOR_PAIREDEND, 
                Double.toString((Integer)configData.get(HS2_MIN_FRAG_LEN_FOR_PAIREDEND)), 
                "1", "NA", logger);
        if(chkPair!=null)
            this.setMinFraglenPairedAlign((Integer)configData.get(HS2_MIN_FRAG_LEN_FOR_PAIREDEND));
        
        chkPair = checkParameter("Integer", HS2_MAX_FRAG_LEN_FOR_PAIREDEND, 
                Double.toString((Integer)configData.get(HS2_MAX_FRAG_LEN_FOR_PAIREDEND)), 
                "1", "NA", logger);
        if(chkPair!=null)
            this.setMaxFraglenPairedAlign((Integer)configData.get(HS2_MAX_FRAG_LEN_FOR_PAIREDEND));
        
        if((Boolean)configData.get(HS2_MATE_ORIENTATION_FR)!=null)
            this.setMateOrientationFR(true);
        if((Boolean)configData.get(HS2_MATE_ORIENTATION_RF)!=null)
            this.setMateOrientationRF(true);
        if((Boolean)configData.get(HS2_MATE_ORIENTATION_FF)!=null)
            this.setMateOrientationFF(true);
        if((Boolean)configData.get(HS2_NO_MIXED_ALN)!=null)
            this.setNoMixedAligns(true);
        if((Boolean)configData.get(HS2_NO_DISCORDANT_ALN)!=null)
            this.setDiscordantAligns(true);
        
        
        logger.info(indent + "passed");
        
    }
    
    
    /**
     * parse HISAT2 parameters defining input settings
     * 
     * @param configData
     * @throws Exception 
     */
    public void parseHisat2OutputConfigurationData(HashMap configData) throws Exception{
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Output parameters");

        if((Boolean)configData.get(HS2_PRINT_WALL_CLOCK)!=null)
            this.setWallClock(true);

        if((Boolean)configData.get(HS2_WRITE_UNPAIRED_READS_NOALN)!=null)
            this.setWriteUnpairedReads(true);
        if((Boolean)configData.get(HS2_WRITE_UNPAIRED_READS_NOALN_GZ)!=null){
            this.setWriteUnpairedReads(true);
            this.setGzipReadFiles(true);
        }
        if((Boolean)configData.get(HS2_WRITE_UNPAIRED_READS_NOALN_BZ2)!=null){
            this.setWriteUnpairedReads(true);
            this.setBz2ReadFiles(true);
        }

        if((Boolean)configData.get(HS2_WRITE_UNPAIRED_READS_ALN)!=null)
            this.setWriteAlignedReads(true);
        if((Boolean)configData.get(HS2_WRITE_UNPAIRED_READS_ALN_GZ)!=null){
            this.setWriteAlignedReads(true);
            this.setGzipReadFiles(true);
        }
        if((Boolean)configData.get(HS2_WRITE_UNPAIRED_READS_ALN_BZ2)!=null){
            this.setWriteAlignedReads(true);
            this.setBz2ReadFiles(true);
        }

        if((Boolean)configData.get(HS2_WRITE_UNPAIRED_READS_DISC_ALN)!=null)
            this.setWriteNonconcordantReads(true);
        if((Boolean)configData.get(HS2_WRITE_UNPAIRED_READS_DISC_ALN_GZ)!=null){
            this.setWriteNonconcordantReads(true);
            this.setGzipReadFiles(true);
        }
        if((Boolean)configData.get(HS2_WRITE_UNPAIRED_READS_DISC_ALN_BZ2)!=null){
            this.setWriteNonconcordantReads(true);
            this.setBz2ReadFiles(true);
        }

        if((Boolean)configData.get(HS2_WRITE_UNPAIRED_READS_CONC_ALN)!=null)
            this.setWriteConcordantReads(true);
        if((Boolean)configData.get(HS2_WRITE_UNPAIRED_READS_CONC_ALN_GZ)!=null){
            this.setWriteConcordantReads(true);
            this.setGzipReadFiles(true);
        }
        if((Boolean)configData.get(HS2_WRITE_UNPAIRED_READS_CONC_ALN_BZ2)!=null){
            this.setWriteConcordantReads(true);
            this.setBz2ReadFiles(true);
        }

        logger.info(indent + "passed");
        
    }
    
    
    /**
     * parse HISAT2 parameters defining input settings
     * 
     * @param configData
     * @throws Exception 
     */
    public void parseHisat2SAMConfigurationData(HashMap configData) throws Exception{
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 SAM parameters");

        if((Boolean)configData.get(HS2_SAM_SUPPRESS_NOALN)!=null)
            this.setSamNoUnaligedReads(true);
        if((Boolean)configData.get(HS2_SAM_SUPPRESS_HEADER_LINES)!=null)
            this.setSamNoHeaderLines(true);
        if((Boolean)configData.get(HS2_SAM_SUPPRESS_SQ_LINES)!=null)
            this.setSamNoSeqHeaderLines(true);

        if(configData.get(HS2_SAM_SET_READ_GROUP_ID)!=null)
            this.setSamReadGroupID((String)configData.get(HS2_SAM_SET_READ_GROUP_ID));
        if(configData.get(HS2_SAM_ADD_READ_GROUP_TEXT)!=null)
            this.setSamReadGroupText((String)configData.get(HS2_SAM_ADD_READ_GROUP_TEXT));

        if(configData.get(HS2_SAM_REM_CHR_TEXT)!=null)
            this.setSamRemoveChrStringFromRef(true);
        if(configData.get(HS2_SAM_ADD_CHR_TEXT)!=null)
            this.setSamAddChrStringFromRef(true);
        if(configData.get(HS2_SAM_OMIT_SEC_AND_QUAL_STRINGS)!=null)
            this.setSamOmitSeqAndQualStrin(true);

        logger.info(indent + "passed");
        
    }
    
    
    /**
     * parse HISAT2 parameters defining input settings
     * 
     * @param configData
     * @throws Exception 
     */
    public void parseHisat2PerformanceConfigurationData(HashMap configData) throws Exception{
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Performance parameters");

        String chkPerf = "";
        chkPerf = checkParameter("Integer", HS2_OVERRIDE_OFFRATE, 
                Integer.toString((Integer)configData.get(HS2_OVERRIDE_OFFRATE)), 
                "1", "NA", logger);
        if(chkPerf!=null)
            this.setIndexOffrate((Integer)configData.get(HS2_OVERRIDE_OFFRATE));

        chkPerf = checkParameter("Integer", HS2_NTHREADS, 
                Integer.toString((Integer)configData.get(HS2_NTHREADS)), 
                "1", "NA", logger);
        if(chkPerf!=null)
            this.setNumOfThreads((Integer)configData.get(HS2_NTHREADS));
        if(configData.get(HS2_REORDER)!=null)
            this.setReorderSAMRecords(true);
        if(configData.get(HS2_MEMORY_MAPPED_IO)!=null)
            this.setUseMemoryMappedIO(true);


        logger.info(indent + "passed");
        
    }
    
    
    /**
     * parse HISAT2 parameters defining input settings
     * 
     * @param configData
     * @throws Exception 
     */
    public void parseHisat2OtherConfigurationData(HashMap configData) throws Exception{
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");

        if(configData.get(HS2_QC_FILTER)!=null)
            this.setQcFiltering(true);
        if(configData.get(HS2_NON_DETERMINISTIC)!=null)
            this.setNonDeterministicAlignment(true);
        
        String chkPerf = "";
        chkPerf = checkParameter("Integer", HS2_SEED, 
                Integer.toString((Integer)configData.get(HS2_SEED)), 
                "0", "NA", logger);
        if(chkPerf!=null)
            this.setSeed((Integer)configData.get(HS2_SEED));

        logger.info(indent + "passed");
        
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
                    Usage
                
                hisat2 [options]* -x <hisat2-idx> 
                 -1 <m1> -2 <m2> 
                 [-S <hit>]

                 */
                this.mapReadsToGenome(sampleData);

                /*
                    write out mapping summary
                */
            } catch (IOException | InterruptedException ex) {
                logger.error("error executing Bowtie Mapping command\n");
                logger.error(cmd);
                logger.error(ex.toString());
            }
        }

    }

    
    
    /**
     * Usage
     * hisat2 [options]* -x <hisat2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <hit>]
     * 
     * @param sampleData
     * @throws IOException
     * @throws InterruptedException 
     */
    private void mapReadsToGenome(SampleDataEntry sampleData) throws IOException, InterruptedException{
                 
        String pathToGenomeIndex = getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + ReferenceDataLocations.ID_REL_BOWTIE_PATH;
                
        ArrayList cmd = new ArrayList<>();
        
        cmd.add(this.getMappingSoftware());
        cmd.add(pathToGenomeIndex);
        
        // add options
        cmd.add(buildOptionsString());
        
        // we assume fastq files 
        cmd.add("-q -1 ");
        cmd.add(sampleData.getFastqFile1());
        if(sampleData.getFastqFile2()!= null)
            cmd.add("-2 " + sampleData.getFastqFile2());
        
        cmd.add("--best");
        cmd.add("-m 2");

        String samGenomeAln = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", SAM_GENALN_EXTENSION);
        cmd.add("--al " + fastqGenomeAln);
        cmd.add("--un " + fastqGenomeUnAln);
        cmd.add("--sam " + samGenomeAln);

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
    
    /**
     * parse the parameters specified in the YAML configuration file to 
     * build the HISAT2 options
     * 
     * because of the number of options, this is broken down into the same
     * subsections specified in the HISAT2 manual
     * 
     * @return 
     */
    private String buildOptionsString(){
        int level = 2;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");
        ArrayList options = new ArrayList<>();
        String optionsString = StringUtils.join(options, " ");
        optionsString = optionsString.concat(this.buildInputOptionsString());
        optionsString = optionsString.concat(this.buildAlignmentOptionsString());
        optionsString = optionsString.concat(this.buildScoringOptionsString());
        optionsString = optionsString.concat(this.buildSpliceOptionsString());
        optionsString = optionsString.concat(this.buildReportingOptionsString());
        optionsString = optionsString.concat(this.buildPairedEndOptionsString());
        optionsString = optionsString.concat(this.buildSAMOptionsString());
        optionsString = optionsString.concat(this.buildPerformanceOptionsString());
        optionsString = optionsString.concat(this.buildOtherOptionsString());
        
        return optionsString;
    }
    
    
    /**
     * build options
     * @return 
     */
    private String buildInputOptionsString(){
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");
        ArrayList options = new ArrayList<>();
        
        String optionsString = "";
        if(this.getSkipFirstNReads()>0)
            options.add("--" + HS2_SKIP_FIRST_NREADS + " "  + Integer.toString(this.getSkipFirstNReads()));
 
        if(this.getAlignFirstNReads()>0)
            options.add("--" + HS2_ALIGN_FIRST_NREADS + " "  + Integer.toString(this.getAlignFirstNReads()));
        
        if(this.getTrimNBasesFrom5Prime()>0)
            options.add("--" + HS2_TRIM_BASES_FROM_5PRIME_END + " "  + Integer.toString(this.getTrimNBasesFrom5Prime()));
        
        if(this.getTrimNBasesFrom3Prime()>0)
            options.add("--" + HS2_TRIM_BASES_FROM_3PRIME_END + " "  + Integer.toString(this.getTrimNBasesFrom3Prime()));

        optionsString = StringUtils.join(options, " ");
        
        
        return optionsString;
        
    }
    
    
    /**
     * build options
     * @return 
     */
    private String buildAlignmentOptionsString(){
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");
        ArrayList options = new ArrayList<>();

        String optionsString = "";
        if(getnCeilingFn().isEmpty()==false & getnCeilingConst()>0 & getnCeilingCoeff()>0 ){
            options.add("--n-ceil");            
            options.add(getnCeilingFn().trim());
            options.add("--" + HS2_N_CEILING_C + " "  + Double.toString(getnCeilingConst()));
            options.add("--" + HS2_N_CEILING_E + " "  + Double.toString(getnCeilingCoeff()));
            optionsString = StringUtils.join(options, " ");        
        }
        return optionsString;
        
    }
    
    
    /**
     * build Scoring options
     * @return 
     */
    private String buildScoringOptionsString(){
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");
        ArrayList options = new ArrayList<>();
        
        options.add("--" + HS2_N_CEILING_C + " "  + Double.toString(getnCeilingConst()));
        
        String optionsString = StringUtils.join(options, " ");
        return optionsString;
        
    }
    
    
    /**
     * build options
     * @return 
     */
    private String buildSpliceOptionsString(){
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");
        ArrayList options = new ArrayList<>();
        
        String optionsString = StringUtils.join(options, " ");
        return optionsString;
        
    }
    
    
    /**
     * build options
     * @return 
     */
    private String buildReportingOptionsString(){
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");
        ArrayList options = new ArrayList<>();
        
        String optionsString = StringUtils.join(options, " ");
        return optionsString;
        
    }
    
    
    /**
     * build options
     * @return 
     */
    private String buildPairedEndOptionsString(){
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");
        ArrayList options = new ArrayList<>();
        
        String optionsString = StringUtils.join(options, " ");
        return optionsString;
        
    }
    
    
    /**
     * build options
     * @return 
     */
    private String buildSAMOptionsString(){
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");
        ArrayList options = new ArrayList<>();
        
        String optionsString = StringUtils.join(options, " ");
        return optionsString;
        
    }
    
    
    /**
     * build options
     * @return 
     */
    private String buildPerformanceOptionsString(){
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");
        ArrayList options = new ArrayList<>();
        
        String optionsString = StringUtils.join(options, " ");
        return optionsString;
        
    }
    
    
    /**
     * build options
     * @return 
     */
    private String buildOtherOptionsString(){
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");
        ArrayList options = new ArrayList<>();
        
        String optionsString = StringUtils.join(options, " ");
        return optionsString;
        
    }
    
    
    @Override
    public void verifyInputData()  throws IOException, NullPointerException{
        logger.info("verify input data");
        
        Validate.notNull((String) this.getMappingSoftware());
        
        if (this.getNumOfThreads() <= 0)
        {
            logger.error("number of threads <" + this.getNumOfThreads() + "> must be positive");    
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

        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + STEP_ID_STRING + ": generate example configuration data");

        HashMap configData = new HashMap();

        configData.put(ID_SOFTWARE, "/usr/local/bin/hisat2");
        configData.put(ID_REF_GENOME, "hsa");

   
        configData.putAll(this.generateExampleInputConfigurationData());
        configData.putAll(this.generateExampleAlignmentConfigurationData());
        configData.putAll(this.generateExampleScoringConfigurationData());
        configData.putAll(this.generateExampleSplicedConfigurationData());
        configData.putAll(this.generateExampleReportingConfigurationData());
        configData.putAll(this.generateExamplePairedEndConfigurationData());
        configData.putAll(this.generateExampleOutputConfigurationData());
        configData.putAll(this.generateExampleSAMConfigurationData());
        configData.putAll(this.generateExamplePerformanceConfigurationData());
        configData.putAll(this.generateExampleOtherConfigurationData());
        

        configData.putAll(configData);
        return configData;
    }

    /**
     * generate sample input configuration data for YAML file
     * @return 
     */
    public HashMap generateExampleInputConfigurationData() {
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");

        HashMap configData = new HashMap();
        configData.put(HS2_SKIP_FIRST_NREADS, "10");
        configData.put(HS2_ALIGN_FIRST_NREADS, "20");
        configData.put(HS2_TRIM_BASES_FROM_5PRIME_END, "30");
        configData.put(HS2_TRIM_BASES_FROM_3PRIME_END, "40");
        return configData;
    }    
    
    
    /**
     * generate sample alignment configuration data for YAML file
     * @return 
     */
    public HashMap generateExampleAlignmentConfigurationData() {
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");

        HashMap configData = new HashMap();
        configData.put(HS2_N_CEILING_F, "10");
        configData.put(HS2_N_CEILING_C, "10");
        configData.put(HS2_N_CEILING_E, "10");
        
        configData.put(HS2_NOFW, "false");
        configData.put(HS2_NOFC, "false");
        configData.put(HS2_IGNORE_QUALITY_VALUES, "false");
        
        return configData;
    }    
    
    
    /**
     * generate sample scoring configuration data for YAML file
     * @return 
     */
    public HashMap generateExampleScoringConfigurationData() {
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");

        HashMap configData = new HashMap();
        configData.put(HS2_PENALTY_MAX_MISMATCH, "21");
        configData.put(HS2_PENALTY_MIN_MISMATCH, "22");
        configData.put(HS2_PENALTY_MAX_SOFTCLIP_MISMATCH, "23");
        configData.put(HS2_PENALTY_MIN_SOFTCLIP_MISMATCH, "24");
        configData.put(HS2_NO_SOFTCLIP, "false");
        configData.put(HS2_PENALTY_AMBIGUOUS_CHARACTER, "25");
        configData.put(HS2_PENALTY_READGAP_OPEN, "26");
        configData.put(HS2_PENALTY_READGAP_EXTEND, "27");
        configData.put(HS2_PENALTY_REFGAP_OPEN, "28");
        configData.put(HS2_PENALTY_REFGAP_EXTEND, "29");
        configData.put(HS2_MIN_ALN_SCORE_FN_C, "L");
        configData.put(HS2_MIN_ALN_SCORE_FN_L, "31.1");
        configData.put(HS2_MIN_ALN_SCORE_FN_S, "31.2");
        return configData;
    }    
    
    
    /**
     * generate sample spliced configuration data for YAML file
     * @return 
     */
    public HashMap generateExampleSplicedConfigurationData() {
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");

        HashMap configData = new HashMap();
        configData.put(HS2_PENALTY_CANONICAL_SPLICE, "40");
        configData.put(HS2_PENALTY_NONCANON_SPLICE, "41");
        configData.put(HS2_PENALTY_LONGINTRON_CAN_SPLICE_F, "A");
        configData.put(HS2_PENALTY_LONGINTRON_CAN_SPLICE_C, "42.1");
        configData.put(HS2_PENALTY_LONGINTRON_CAN_SPLICE_E, "42.2");
        configData.put(HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_F, "B");
        configData.put(HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_C, "43.1");
        configData.put(HS2_PENALTY_LONGINTRON_NONCAN_SPLICE_E, "43.2");
        configData.put(HS2_MIN_INTRON_LEN, "44");
        configData.put(HS2_MAX_INTRON_LEN, "45");
        configData.put(HS2_NO_SPLICED_ALIGNMENT, "true");
        configData.put(HS2_TRANSCRIPTOME_MAPPING_ONLY, "false");
        configData.put(HS2_DOWNSTREAM_TRANSCRIPTOME_ASSEMBLY, "true");
        configData.put(HS2_DTA_CUFFLINKS, "false");
        configData.put(HS2_KNOWN_SPLICES_SITES_INFILE, "KnownSpliceIn.txt");
        configData.put(HS2_NOVEL_SPLICES_SITES_INFILE, "NovelSpliceIn.txt");
        configData.put(HS2_NOVEL_SPLICES_SITES_OUTFILE, "NovelSpliceOut.txt");
        configData.put(HS2_RNA_STRANDEDNESS, "RF");
        return configData;
    }    
    
    
    /**
     * generate sample reporting configuration data for YAML file
     * @return 
     */
    public HashMap generateExampleReportingConfigurationData() {
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");

        HashMap configData = new HashMap();
        configData.put(HS2_NUM_DISTINCT_PRIMARY_ALNS, "50");
        configData.put(HS2_MAX_NUM_SEEDS, "51");
        configData.put(HS2_REPORT_SECONDARY_ALNS, "false");
        return configData;
    }    
    
    
    /**
     * generate sample paired end configuration data for YAML file
     * @return 
     */
    public HashMap generateExamplePairedEndConfigurationData() {
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");

        HashMap configData = new HashMap();
        configData.put(HS2_MIN_FRAG_LEN_FOR_PAIREDEND, "60");
        configData.put(HS2_MAX_FRAG_LEN_FOR_PAIREDEND, "61");
        configData.put(HS2_MATE_ORIENTATION_FR, "true");
        configData.put(HS2_MATE_ORIENTATION_RF, "false");
        configData.put(HS2_MATE_ORIENTATION_FF, "false");
        configData.put(HS2_NO_MIXED_ALN, "false");
        configData.put(HS2_NO_DISCORDANT_ALN, "false");
        return configData;
    }    
    
    
    /**
     * generate sample output configuration data for YAML file
     * @return 
     */
    public HashMap generateExampleOutputConfigurationData() {
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");

        HashMap configData = new HashMap();
        configData.put(HS2_NO_DISCORDANT_ALN, "false");
        configData.put(HS2_WRITE_UNPAIRED_READS_NOALN, "true");
        configData.put(HS2_WRITE_UNPAIRED_READS_NOALN_GZ, "false");
        configData.put(HS2_WRITE_UNPAIRED_READS_NOALN_BZ2, "true");
        configData.put(HS2_WRITE_UNPAIRED_READS_ALN, "true");
        configData.put(HS2_WRITE_UNPAIRED_READS_ALN_GZ, "false");
        configData.put(HS2_WRITE_UNPAIRED_READS_ALN_BZ2, "true");
        configData.put(HS2_WRITE_UNPAIRED_READS_DISC_ALN, "true");
        configData.put(HS2_WRITE_UNPAIRED_READS_DISC_ALN_GZ, "false");
        configData.put(HS2_WRITE_UNPAIRED_READS_DISC_ALN_BZ2, "true");
        configData.put(HS2_WRITE_UNPAIRED_READS_CONC_ALN, "true");
        configData.put(HS2_WRITE_UNPAIRED_READS_CONC_ALN_GZ, "false");
        configData.put(HS2_WRITE_UNPAIRED_READS_CONC_ALN_BZ2, "true");
        return configData;
    }    
    
    
    /**
     * generate sample SAM configuration data for YAML file
     * @return 
     */
    public HashMap generateExampleSAMConfigurationData() {
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");

        HashMap configData = new HashMap();
        configData.put(HS2_SAM_SUPPRESS_NOALN, "false");
        configData.put(HS2_SAM_SUPPRESS_HEADER_LINES, "false");
        configData.put(HS2_SAM_SUPPRESS_SQ_LINES, "true");
        configData.put(HS2_SAM_SET_READ_GROUP_ID, "ReadGroupID");
        configData.put(HS2_SAM_ADD_READ_GROUP_TEXT, "ReadGroupText");
        configData.put(HS2_SAM_REM_CHR_TEXT, "false");
        configData.put(HS2_SAM_ADD_CHR_TEXT, "true");
        configData.put(HS2_SAM_OMIT_SEC_AND_QUAL_STRINGS, "false");
        return configData;
    }    
    
    
    /**
     * generate sample performance configuration data for YAML file
     * @return 
     */
    public HashMap generateExamplePerformanceConfigurationData() {
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");

        HashMap configData = new HashMap();
        configData.put(HS2_OVERRIDE_OFFRATE, "false");
        configData.put(HS2_NTHREADS, "4");
        configData.put(HS2_REORDER, "false");
        configData.put(HS2_MEMORY_MAPPED_IO, "false");
        return configData;
    }    
    
    
    /**
     * generate sample other configuration data for YAML file
     * @return 
     */
    public HashMap generateExampleOtherConfigurationData() {
        int level = 3;
        String indent = StringUtils.repeat("-", level*2);
        logger.info(indent + "parsing HISAT2 Other parameters");

        HashMap configData = new HashMap();
        configData.put(HS2_QC_FILTER, "false");
        configData.put(HS2_NON_DETERMINISTIC, "false");
        configData.put(HS2_SEED, "70");
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
     * @return the nCeilingFn
     */
    public String getnCeilingFn() {
        return nCeilingFn;
    }

    /**
     * @param nCeilingFn the nCeilingFn to set
     */
    public void setnCeilingFn(String nCeilingFn) {
        this.nCeilingFn = nCeilingFn;
    }

    /**
     * @return the nCeilingConst
     */
    public double getnCeilingConst() {
        return nCeilingConst;
    }

    /**
     * @param nCeilingConst the nCeilingConst to set
     */
    public void setnCeilingConst(double nCeilingConst) {
        this.nCeilingConst = nCeilingConst;
    }

    /**
     * @return the nCeilingCoeff
     */
    public double getnCeilingCoeff() {
        return nCeilingCoeff;
    }

    /**
     * @param nCeilingCoeff the nCeilingCoeff to set
     */
    public void setnCeilingCoeff(double nCeilingCoeff) {
        this.nCeilingCoeff = nCeilingCoeff;
    }

    /**
     * @return the nofw
     */
    public Boolean getNofw() {
        return nofw;
    }

    /**
     * @param nofw the nofw to set
     */
    public void setNofw(Boolean nofw) {
        this.nofw = nofw;
    }

    /**
     * @return the nofc
     */
    public Boolean getNofc() {
        return nofc;
    }

    /**
     * @param nofc the nofc to set
     */
    public void setNofc(Boolean nofc) {
        this.nofc = nofc;
    }

    /**
     * @return the ignoreQuals
     */
    public Boolean getIgnoreQuals() {
        return ignoreQuals;
    }

    /**
     * @param ignoreQuals the ignoreQuals to set
     */
    public void setIgnoreQuals(Boolean ignoreQuals) {
        this.ignoreQuals = ignoreQuals;
    }

    /**
     * @return the skipFirstNReads
     */
    public int getSkipFirstNReads() {
        return skipFirstNReads;
    }

    /**
     * @param skipFirstNReads the skipFirstNReads to set
     */
    public void setSkipFirstNReads(int skipFirstNReads) {
        this.skipFirstNReads = skipFirstNReads;
    }

    /**
     * @return the alignFirstNReads
     */
    public int getAlignFirstNReads() {
        return alignFirstNReads;
    }

    /**
     * @param alignFirstNReads the alignFirstNReads to set
     */
    public void setAlignFirstNReads(int alignFirstNReads) {
        this.alignFirstNReads = alignFirstNReads;
    }

    /**
     * @return the trimNBasesFrom5Prime
     */
    public int getTrimNBasesFrom5Prime() {
        return trimNBasesFrom5Prime;
    }

    /**
     * @param trimNBasesFrom5Prime the trimNBasesFrom5Prime to set
     */
    public void setTrimNBasesFrom5Prime(int trimNBasesFrom5Prime) {
        this.trimNBasesFrom5Prime = trimNBasesFrom5Prime;
    }

    /**
     * @return the trimNBasesFrom3Prime
     */
    public int getTrimNBasesFrom3Prime() {
        return trimNBasesFrom3Prime;
    }

    /**
     * @param trimNBasesFrom3Prime the trimNBasesFrom3Prime to set
     */
    public void setTrimNBasesFrom3Prime(int trimNBasesFrom3Prime) {
        this.trimNBasesFrom3Prime = trimNBasesFrom3Prime;
    }

    /**
     * @return the numOfThreads
     */
    public int getNumOfThreads() {
        return numOfThreads;
    }

    /**
     * @param numOfThreads the numOfThreads to set
     */
    public void setNumOfThreads(int numOfThreads) {
        this.numOfThreads = numOfThreads;
    }

    /**
     * @return the maxMismatchPenalty
     */
    public int getMaxMismatchPenalty() {
        return maxMismatchPenalty;
    }

    /**
     * @param maxMismatchPenalty the maxMismatchPenalty to set
     */
    public void setMaxMismatchPenalty(int maxMismatchPenalty) {
        this.maxMismatchPenalty = maxMismatchPenalty;
    }

    /**
     * @return the minMismatchPenalty
     */
    public int getMinMismatchPenalty() {
        return minMismatchPenalty;
    }

    /**
     * @param minMismatchPenalty the minMismatchPenalty to set
     */
    public void setMinMismatchPenalty(int minMismatchPenalty) {
        this.minMismatchPenalty = minMismatchPenalty;
    }

    /**
     * @return the maxsoftClipPenalty
     */
    public int getMaxsoftClipPenalty() {
        return maxsoftClipPenalty;
    }

    /**
     * @param maxsoftClipPenalty the maxsoftClipPenalty to set
     */
    public void setMaxsoftClipPenalty(int maxsoftClipPenalty) {
        this.maxsoftClipPenalty = maxsoftClipPenalty;
    }

    /**
     * @return the minsoftClipPenalty
     */
    public int getMinsoftClipPenalty() {
        return minsoftClipPenalty;
    }

    /**
     * @param minsoftClipPenalty the minsoftClipPenalty to set
     */
    public void setMinsoftClipPenalty(int minsoftClipPenalty) {
        this.minsoftClipPenalty = minsoftClipPenalty;
    }

    /**
     * @return the allowSoftClip
     */
    public Boolean getAllowSoftClip() {
        return allowSoftClip;
    }

    /**
     * @param allowSoftClip the allowSoftClip to set
     */
    public void setAllowSoftClip(Boolean allowSoftClip) {
        this.allowSoftClip = allowSoftClip;
    }

    /**
     * @return the penaltyAmbiguousChar
     */
    public int getPenaltyAmbiguousChar() {
        return penaltyAmbiguousChar;
    }

    /**
     * @param penaltyAmbiguousChar the penaltyAmbiguousChar to set
     */
    public void setPenaltyAmbiguousChar(int penaltyAmbiguousChar) {
        this.penaltyAmbiguousChar = penaltyAmbiguousChar;
    }

    /**
     * @return the readGapOpenPenalty
     */
    public int getReadGapOpenPenalty() {
        return readGapOpenPenalty;
    }

    /**
     * @param readGapOpenPenalty the readGapOpenPenalty to set
     */
    public void setReadGapOpenPenalty(int readGapOpenPenalty) {
        this.readGapOpenPenalty = readGapOpenPenalty;
    }

    /**
     * @return the readGapExtendPenalty
     */
    public int getReadGapExtendPenalty() {
        return readGapExtendPenalty;
    }

    /**
     * @param readGapExtendPenalty the readGapExtendPenalty to set
     */
    public void setReadGapExtendPenalty(int readGapExtendPenalty) {
        this.readGapExtendPenalty = readGapExtendPenalty;
    }

    /**
     * @return the refGapOpenPenalty
     */
    public int getRefGapOpenPenalty() {
        return refGapOpenPenalty;
    }

    /**
     * @param refGapOpenPenalty the refGapOpenPenalty to set
     */
    public void setRefGapOpenPenalty(int refGapOpenPenalty) {
        this.refGapOpenPenalty = refGapOpenPenalty;
    }

    /**
     * @return the refGapExtendPenalty
     */
    public int getRefGapExtendPenalty() {
        return refGapExtendPenalty;
    }

    /**
     * @param refGapExtendPenalty the refGapExtendPenalty to set
     */
    public void setRefGapExtendPenalty(int refGapExtendPenalty) {
        this.refGapExtendPenalty = refGapExtendPenalty;
    }

    /**
     * @return the scoreMinFun
     */
    public String getScoreMinFun() {
        return scoreMinFun;
    }

    /**
     * @param scoreMinFun the scoreMinFun to set
     */
    public void setScoreMinFun(String scoreMinFun) {
        this.scoreMinFun = scoreMinFun;
    }

    /**
     * @return the scoreMinConst
     */
    public double getScoreMinConst() {
        return scoreMinConst;
    }

    /**
     * @param scoreMinConst the scoreMinConst to set
     */
    public void setScoreMinConst(double scoreMinConst) {
        this.scoreMinConst = scoreMinConst;
    }

    /**
     * @return the scoreMinCoeff
     */
    public double getScoreMinCoeff() {
        return scoreMinCoeff;
    }

    /**
     * @param scoreMinCoeff the scoreMinCoeff to set
     */
    public void setScoreMinCoeff(double scoreMinCoeff) {
        this.scoreMinCoeff = scoreMinCoeff;
    }

    /**
     * @return the penaltyCanSplice
     */
    public int getPenaltyCanSplice() {
        return penaltyCanSplice;
    }

    /**
     * @param penaltyCanSplice the penaltyCanSplice to set
     */
    public void setPenaltyCanSplice(int penaltyCanSplice) {
        this.penaltyCanSplice = penaltyCanSplice;
    }

    /**
     * @return the penaltyNonCanSplice
     */
    public int getPenaltyNonCanSplice() {
        return penaltyNonCanSplice;
    }

    /**
     * @param penaltyNonCanSplice the penaltyNonCanSplice to set
     */
    public void setPenaltyNonCanSplice(int penaltyNonCanSplice) {
        this.penaltyNonCanSplice = penaltyNonCanSplice;
    }

    /**
     * @return the penaltyLongCanSpiceFun
     */
    public String getPenaltyLongCanSpiceFun() {
        return penaltyLongCanSpiceFun;
    }

    /**
     * @param penaltyLongCanSpiceFun the penaltyLongCanSpiceFun to set
     */
    public void setPenaltyLongCanSpiceFun(String penaltyLongCanSpiceFun) {
        this.penaltyLongCanSpiceFun = penaltyLongCanSpiceFun;
    }

    /**
     * @return the penaltyLongCanSpiceConst
     */
    public double getPenaltyLongCanSpiceConst() {
        return penaltyLongCanSpiceConst;
    }

    /**
     * @param penaltyLongCanSpiceConst the penaltyLongCanSpiceConst to set
     */
    public void setPenaltyLongCanSpiceConst(double penaltyLongCanSpiceConst) {
        this.penaltyLongCanSpiceConst = penaltyLongCanSpiceConst;
    }

    /**
     * @return the penaltyLongCanSpiceCoeff
     */
    public double getPenaltyLongCanSpiceCoeff() {
        return penaltyLongCanSpiceCoeff;
    }

    /**
     * @param penaltyLongCanSpiceCoeff the penaltyLongCanSpiceCoeff to set
     */
    public void setPenaltyLongCanSpiceCoeff(double penaltyLongCanSpiceCoeff) {
        this.penaltyLongCanSpiceCoeff = penaltyLongCanSpiceCoeff;
    }

    /**
     * @return the penaltyLongNonCanSpiceFun
     */
    public String getPenaltyLongNonCanSpiceFun() {
        return penaltyLongNonCanSpiceFun;
    }

    /**
     * @param penaltyLongNonCanSpiceFun the penaltyLongNonCanSpiceFun to set
     */
    public void setPenaltyLongNonCanSpiceFun(String penaltyLongNonCanSpiceFun) {
        this.penaltyLongNonCanSpiceFun = penaltyLongNonCanSpiceFun;
    }

    /**
     * @return the penaltyLongNonCanSpiceConst
     */
    public double getPenaltyLongNonCanSpiceConst() {
        return penaltyLongNonCanSpiceConst;
    }

    /**
     * @param penaltyLongNonCanSpiceConst the penaltyLongNonCanSpiceConst to set
     */
    public void setPenaltyLongNonCanSpiceConst(double penaltyLongNonCanSpiceConst) {
        this.penaltyLongNonCanSpiceConst = penaltyLongNonCanSpiceConst;
    }

    /**
     * @return the penaltyLongNonCanSpiceCoeff
     */
    public double getPenaltyLongNonCanSpiceCoeff() {
        return penaltyLongNonCanSpiceCoeff;
    }

    /**
     * @param penaltyLongNonCanSpiceCoeff the penaltyLongNonCanSpiceCoeff to set
     */
    public void setPenaltyLongNonCanSpiceCoeff(double penaltyLongNonCanSpiceCoeff) {
        this.penaltyLongNonCanSpiceCoeff = penaltyLongNonCanSpiceCoeff;
    }

    /**
     * @return the minIntronLength
     */
    public int getMinIntronLength() {
        return minIntronLength;
    }

    /**
     * @param minIntronLength the minIntronLength to set
     */
    public void setMinIntronLength(int minIntronLength) {
        this.minIntronLength = minIntronLength;
    }

    /**
     * @return the maxIntronLength
     */
    public int getMaxIntronLength() {
        return maxIntronLength;
    }

    /**
     * @param maxIntronLength the maxIntronLength to set
     */
    public void setMaxIntronLength(int maxIntronLength) {
        this.maxIntronLength = maxIntronLength;
    }

    /**
     * @return the knownSpliceSitesFileIn
     */
    public String getKnownSpliceSitesFileIn() {
        return knownSpliceSitesFileIn;
    }

    /**
     * @param knownSpliceSitesFileIn the knownSpliceSitesFileIn to set
     */
    public void setKnownSpliceSitesFileIn(String knownSpliceSitesFileIn) {
        this.knownSpliceSitesFileIn = knownSpliceSitesFileIn;
    }

    /**
     * @return the novelSpliceSitesFileIn
     */
    public String getNovelSpliceSitesFileIn() {
        return novelSpliceSitesFileIn;
    }

    /**
     * @param novelSpliceSitesFileIn the novelSpliceSitesFileIn to set
     */
    public void setNovelSpliceSitesFileIn(String novelSpliceSitesFileIn) {
        this.novelSpliceSitesFileIn = novelSpliceSitesFileIn;
    }

    /**
     * @return the novelSpliceSitesFileOut
     */
    public String getNovelSpliceSitesFileOut() {
        return novelSpliceSitesFileOut;
    }

    /**
     * @param novelSpliceSitesFileOut the novelSpliceSitesFileOut to set
     */
    public void setNovelSpliceSitesFileOut(String novelSpliceSitesFileOut) {
        this.novelSpliceSitesFileOut = novelSpliceSitesFileOut;
    }

    /**
     * @return the disableSplicedAligned
     */
    public Boolean getDisableSplicedAligned() {
        return disableSplicedAligned;
    }

    /**
     * @param disableSplicedAligned the disableSplicedAligned to set
     */
    public void setDisableSplicedAligned(Boolean disableSplicedAligned) {
        this.disableSplicedAligned = disableSplicedAligned;
    }

    /**
     * @return the rnaStrandness
     */
    public String getRnaStrandness() {
        return rnaStrandness;
    }

    /**
     * @param rnaStrandness the rnaStrandness to set
     */
    public void setRnaStrandness(String rnaStrandness) {
        this.rnaStrandness = rnaStrandness;
    }

    /**
     * @return the transcriptomeAnalysisOnly
     */
    public Boolean getTranscriptomeAnalysisOnly() {
        return transcriptomeAnalysisOnly;
    }

    /**
     * @param transcriptomeAnalysisOnly the transcriptomeAnalysisOnly to set
     */
    public void setTranscriptomeAnalysisOnly(Boolean transcriptomeAnalysisOnly) {
        this.transcriptomeAnalysisOnly = transcriptomeAnalysisOnly;
    }

    /**
     * @return the dta
     */
    public Boolean getDta() {
        return dta;
    }

    /**
     * @param dta the dta to set
     */
    public void setDta(Boolean dta) {
        this.dta = dta;
    }

    /**
     * @return the dtrCufflinks
     */
    public Boolean getDtaCufflinks() {
        return dtaCufflinks;
    }

    /**
     * @param dtrCufflinks the dtrCufflinks to set
     */
    public void setDtaCufflinks(Boolean dtrCufflinks) {
        this.dtaCufflinks = dtrCufflinks;
    }

    /**
     * @return the disableTemplateAdjustment
     */
    public Boolean getDisableTemplateAdjustment() {
        return disableTemplateAdjustment;
    }

    /**
     * @param disableTemplateAdjustment the disableTemplateAdjustment to set
     */
    public void setDisableTemplateAdjustment(Boolean disableTemplateAdjustment) {
        this.disableTemplateAdjustment = disableTemplateAdjustment;
    }

    /**
     * @return the minFraglenPairedAlign
     */
    public int getMinFraglenPairedAlign() {
        return minFraglenPairedAlign;
    }

    /**
     * @param minFraglenPairedAlign the minFraglenPairedAlign to set
     */
    public void setMinFraglenPairedAlign(int minFraglenPairedAlign) {
        this.minFraglenPairedAlign = minFraglenPairedAlign;
    }

    /**
     * @return the maxFraglenPairedAlign
     */
    public int getMaxFraglenPairedAlign() {
        return maxFraglenPairedAlign;
    }

    /**
     * @param maxFraglenPairedAlign the maxFraglenPairedAlign to set
     */
    public void setMaxFraglenPairedAlign(int maxFraglenPairedAlign) {
        this.maxFraglenPairedAlign = maxFraglenPairedAlign;
    }

    /**
     * @return the noMixedAligns
     */
    public Boolean getNoMixedAligns() {
        return noMixedAligns;
    }

    /**
     * @param noMixedAligns the noMixedAligns to set
     */
    public void setNoMixedAligns(Boolean noMixedAligns) {
        this.noMixedAligns = noMixedAligns;
    }

    /**
     * @return the discordantAligns
     */
    public Boolean getDiscordantAligns() {
        return discordantAligns;
    }

    /**
     * @param discordantAligns the discordantAligns to set
     */
    public void setDiscordantAligns(Boolean discordantAligns) {
        this.discordantAligns = discordantAligns;
    }

    /**
     * @return the numOfDistinctPrimaryAligns
     */
    public int getNumOfDistinctPrimaryAligns() {
        return numOfDistinctPrimaryAligns;
    }

    /**
     * @param numOfDistinctPrimaryAligns the numOfDistinctPrimaryAligns to set
     */
    public void setNumOfDistinctPrimaryAligns(int numOfDistinctPrimaryAligns) {
        this.numOfDistinctPrimaryAligns = numOfDistinctPrimaryAligns;
    }

    /**
     * @return the maxNumSeeds
     */
    public int getMaxNumSeeds() {
        return maxNumSeeds;
    }

    /**
     * @param maxNumSeeds the maxNumSeeds to set
     */
    public void setMaxNumSeeds(int maxNumSeeds) {
        this.maxNumSeeds = maxNumSeeds;
    }

    /**
     * @return the reportSecondaryAligns
     */
    public Boolean getReportSecondaryAligns() {
        return reportSecondaryAligns;
    }

    /**
     * @param reportSecondaryAligns the reportSecondaryAligns to set
     */
    public void setReportSecondaryAligns(Boolean reportSecondaryAligns) {
        this.reportSecondaryAligns = reportSecondaryAligns;
    }

    /**
     * @return the mateOrientationFR
     */
    public Boolean getMateOrientationFR() {
        return mateOrientationFR;
    }

    /**
     * @param mateOrientationFR the mateOrientationFR to set
     */
    public void setMateOrientationFR(Boolean mateOrientationFR) {
        this.mateOrientationFR = mateOrientationFR;
    }

    /**
     * @return the mateOrientationRF
     */
    public Boolean getMateOrientationRF() {
        return mateOrientationRF;
    }

    /**
     * @param mateOrientationRF the mateOrientationRF to set
     */
    public void setMateOrientationRF(Boolean mateOrientationRF) {
        this.mateOrientationRF = mateOrientationRF;
    }

    /**
     * @return the mateOrientationFFF
     */
    public Boolean getMateOrientationFF() {
        return mateOrientationFF;
    }

    /**
     * @param mateOrientationFFF the mateOrientationFFF to set
     */
    public void setMateOrientationFF(Boolean mateOrientationFFF) {
        this.mateOrientationFF = mateOrientationFFF;
    }

    /**
     * @return the quietMode
     */
    public Boolean getQuietMode() {
        return quietMode;
    }

    /**
     * @param quietMode the quietMode to set
     */
    public void setQuietMode(Boolean quietMode) {
        this.quietMode = quietMode;
    }

    /**
     * @return the writeAlignmentSummary
     */
    public Boolean getWriteAlignmentSummary() {
        return writeAlignmentSummary;
    }

    /**
     * @param writeAlignmentSummary the writeAlignmentSummary to set
     */
    public void setWriteAlignmentSummary(Boolean writeAlignmentSummary) {
        this.writeAlignmentSummary = writeAlignmentSummary;
    }

    /**
     * @return the writeAlignmentNewStyle
     */
    public Boolean getWriteAlignmentNewStyle() {
        return writeAlignmentNewStyle;
    }

    /**
     * @param writeAlignmentNewStyle the writeAlignmentNewStyle to set
     */
    public void setWriteAlignmentNewStyle(Boolean writeAlignmentNewStyle) {
        this.writeAlignmentNewStyle = writeAlignmentNewStyle;
    }

    /**
     * @return the writeHisatMetrics
     */
    public Boolean getWriteHisatMetrics() {
        return writeHisatMetrics;
    }

    /**
     * @param writeHisatMetrics the writeHisatMetrics to set
     */
    public void setWriteHisatMetrics(Boolean writeHisatMetrics) {
        this.writeHisatMetrics = writeHisatMetrics;
    }

    /**
     * @return the writeHisatMetricsStderr
     */
    public Boolean getWriteHisatMetricsStderr() {
        return writeHisatMetricsStderr;
    }

    /**
     * @param writeHisatMetricsStderr the writeHisatMetricsStderr to set
     */
    public void setWriteHisatMetricsStderr(Boolean writeHisatMetricsStderr) {
        this.writeHisatMetricsStderr = writeHisatMetricsStderr;
    }

    /**
     * @return the metWriteInterval
     */
    public int getMetWriteInterval() {
        return metWriteInterval;
    }

    /**
     * @param metWriteInterval the metWriteInterval to set
     */
    public void setMetWriteInterval(int metWriteInterval) {
        this.metWriteInterval = metWriteInterval;
    }

    /**
     * @return the samNoUnaligedReads
     */
    public Boolean getSamNoUnaligedReads() {
        return samNoUnaligedReads;
    }

    /**
     * @param samNoUnaligedReads the samNoUnaligedReads to set
     */
    public void setSamNoUnaligedReads(Boolean samNoUnaligedReads) {
        this.samNoUnaligedReads = samNoUnaligedReads;
    }

    /**
     * @return the samNoHeaderLines
     */
    public Boolean getSamNoHeaderLines() {
        return samNoHeaderLines;
    }

    /**
     * @param samNoHeaderLines the samNoHeaderLines to set
     */
    public void setSamNoHeaderLines(Boolean samNoHeaderLines) {
        this.samNoHeaderLines = samNoHeaderLines;
    }

    /**
     * @return the samNoSeqHeaderLines
     */
    public Boolean getSamNoSeqHeaderLines() {
        return samNoSeqHeaderLines;
    }

    /**
     * @param samNoSeqHeaderLines the samNoSeqHeaderLines to set
     */
    public void setSamNoSeqHeaderLines(Boolean samNoSeqHeaderLines) {
        this.samNoSeqHeaderLines = samNoSeqHeaderLines;
    }

    /**
     * @return the samReadGroupID
     */
    public String getSamReadGroupID() {
        return samReadGroupID;
    }

    /**
     * @param samReadGroupID the samReadGroupID to set
     */
    public void setSamReadGroupID(String samReadGroupID) {
        this.samReadGroupID = samReadGroupID;
    }

    /**
     * @return the samReadGroupText
     */
    public String getSamReadGroupText() {
        return samReadGroupText;
    }

    /**
     * @param samReadGroupText the samReadGroupText to set
     */
    public void setSamReadGroupText(String samReadGroupText) {
        this.samReadGroupText = samReadGroupText;
    }

    /**
     * @return the samRemoveChrStringFromRef
     */
    public Boolean getSamRemoveChrStringFromRef() {
        return samRemoveChrStringFromRef;
    }

    /**
     * @param samRemoveChrStringFromRef the samRemoveChrStringFromRef to set
     */
    public void setSamRemoveChrStringFromRef(Boolean samRemoveChrStringFromRef) {
        this.samRemoveChrStringFromRef = samRemoveChrStringFromRef;
    }

    /**
     * @return the samAddChrStringFromRef
     */
    public Boolean getSamAddChrStringFromRef() {
        return samAddChrStringFromRef;
    }

    /**
     * @param samAddChrStringFromRef the samAddChrStringFromRef to set
     */
    public void setSamAddChrStringFromRef(Boolean samAddChrStringFromRef) {
        this.samAddChrStringFromRef = samAddChrStringFromRef;
    }

    /**
     * @return the samOmitSeqAndQualStrin
     */
    public Boolean getSamOmitSeqAndQualStrin() {
        return samOmitSeqAndQualStrin;
    }

    /**
     * @param samOmitSeqAndQualStrin the samOmitSeqAndQualStrin to set
     */
    public void setSamOmitSeqAndQualStrin(Boolean samOmitSeqAndQualStrin) {
        this.samOmitSeqAndQualStrin = samOmitSeqAndQualStrin;
    }

    /**
     * @return the wallClock
     */
    public Boolean getWallClock() {
        return wallClock;
    }

    /**
     * @param wallClock the wallClock to set
     */
    public void setWallClock(Boolean wallClock) {
        this.wallClock = wallClock;
    }

    /**
     * @return the writeUnpairedReads
     */
    public Boolean getWriteUnpairedReads() {
        return writeUnpairedReads;
    }

    /**
     * @param writeUnpairedReads the writeUnpairedReads to set
     */
    public void setWriteUnpairedReads(Boolean writeUnpairedReads) {
        this.writeUnpairedReads = writeUnpairedReads;
    }

    /**
     * @return the writeAlignedReads
     */
    public Boolean getWriteAlignedReads() {
        return writeAlignedReads;
    }

    /**
     * @param writeAlignedReads the writeAlignedReads to set
     */
    public void setWriteAlignedReads(Boolean writeAlignedReads) {
        this.writeAlignedReads = writeAlignedReads;
    }

    /**
     * @return the writeConcordantReads
     */
    public Boolean getWriteConcordantReads() {
        return writeConcordantReads;
    }

    /**
     * @param writeConcordantReads the writeConcordantReads to set
     */
    public void setWriteConcordantReads(Boolean writeConcordantReads) {
        this.writeConcordantReads = writeConcordantReads;
    }

    /**
     * @return the writeNonconcordantReads
     */
    public Boolean getWriteNonconcordantReads() {
        return writeNonconcordantReads;
    }

    /**
     * @param writeNonconcordantReads the writeNonconcordantReads to set
     */
    public void setWriteNonconcordantReads(Boolean writeNonconcordantReads) {
        this.writeNonconcordantReads = writeNonconcordantReads;
    }

    /**
     * @return the gzipReadFiles
     */
    public Boolean getGzipReadFiles() {
        return gzipReadFiles;
    }

    /**
     * @param gzipReadFiles the gzipReadFiles to set
     */
    public void setGzipReadFiles(Boolean gzipReadFiles) {
        this.gzipReadFiles = gzipReadFiles;
    }

    /**
     * @return the bz2ReadFiles
     */
    public Boolean getBz2ReadFiles() {
        return bz2ReadFiles;
    }

    /**
     * @param bz2ReadFiles the bz2ReadFiles to set
     */
    public void setBz2ReadFiles(Boolean bz2ReadFiles) {
        this.bz2ReadFiles = bz2ReadFiles;
    }

    /**
     * @return the indexOffrate
     */
    public int getIndexOffrate() {
        return indexOffrate;
    }

    /**
     * @param indexOffrate the indexOffrate to set
     */
    public void setIndexOffrate(int indexOffrate) {
        this.indexOffrate = indexOffrate;
    }

    /**
     * @return the reorderSAMRecords
     */
    public Boolean getReorderSAMRecords() {
        return reorderSAMRecords;
    }

    /**
     * @param reorderSAMRecords the reorderSAMRecords to set
     */
    public void setReorderSAMRecords(Boolean reorderSAMRecords) {
        this.reorderSAMRecords = reorderSAMRecords;
    }

    /**
     * @return the useMemoryMappedIO
     */
    public Boolean getUseMemoryMappedIO() {
        return useMemoryMappedIO;
    }

    /**
     * @param useMemoryMappedIO the useMemoryMappedIO to set
     */
    public void setUseMemoryMappedIO(Boolean useMemoryMappedIO) {
        this.useMemoryMappedIO = useMemoryMappedIO;
    }

        /**
     * @return the qcFiltering
     */
    public Boolean getQcFiltering() {
        return qcFiltering;
    }

    /**
     * @param qcFiltering the qcFiltering to set
     */
    public void setQcFiltering(Boolean qcFiltering) {
        this.qcFiltering = qcFiltering;
    }

    /**
     * @return the seed
     */
    public int getSeed() {
        return seed;
    }

    /**
     * @param seed the seed to set
     */
    public void setSeed(int seed) {
        this.seed = seed;
    }

    /**
     * @return the nonDeterministicAlignment
     */
    public Boolean getNonDeterministicAlignment() {
        return nonDeterministicAlignment;
    }

    /**
     * @param nonDeterministicAlignment the nonDeterministicAlignment to set
     */
    public void setNonDeterministicAlignment(Boolean nonDeterministicAlignment) {
        this.nonDeterministicAlignment = nonDeterministicAlignment;
    }

}
