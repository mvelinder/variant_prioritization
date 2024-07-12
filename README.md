# variant_prioritization
purpose: a curated list of resources and tools (mostly) for variant prioritization and interpretation in clinical and research settings, always welcoming additions and changes

## web tools / data sets

### aggregate tools
- [Varsome](https://varsome.com/) - aggregate gene and variant information 
- [MARRVEL](http://marrvel.org/) - model organism focused gene and variant information
- [GenCC](https://thegencc.org) - multiple gene:disease curation efforts in one place
- [Franklin](https://franklin.genoox.com/) - variant interpretation and information
- [GeneCards](http://www.genecards.org/) - aggregate gene information
- [PhenCards](https://phencards.org/) - phenotype catalog and searching
- [HumanMine](http://www.humanmine.org/) - aggregate gene information
- [MonarchInitiative](https://monarchinitiative.org/) - gene, disease, biological function information, including gene ortholog phenotypes
- [iobio](http://iobio.io/) - real time genomic analysis tools for QC and variant prioritization
- [UCSCbrowser](https://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=724544049_pngh3ffiA9LYDPiWojHaNAcDu3CA) - supremely customizable genome browser
- [WashUStLEpigenome](http://epigenomegateway.wustl.edu/browser/) - epigenome browser and visualization

### manually curated knowledge bases
- [OMIM](http://omim.org/) - Online Mendelian Inheritance in Man, gene:disease association, gene function, variant molecular genetics
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) - clinical variant interpretations
- [ClinGen](https://www.clinicalgenome.org/) - expert curated gene:disease associations
  - [CNV-Interpretation](http://cnvcalc.clinicalgenome.org/cnvcalc/) - CNV pathogenicity calculator
- [HGMD](http://www.hgmd.cf.ac.uk/ac/index.php) - curated gene and variant disease associations based on literature, with public and paid version
- [MedGen](https://www.ncbi.nlm.nih.gov/medgen/) - clinically focused gene and disease information including literature review
- [SimpleClinVar](http://simple-clinvar.broadinstitute.org/) - simplied summaries and views of ClinVar data

### biological impact predictions
- [MaveDB](https://www.mavedb.org/) - database of Multiplexed Assays of Variant Effect
- [EVE](https://evemodel.org/) - evolutionary model of missense variant effects
- [REVEL](https://sites.google.com/site/revelgenomics/) - missense variant pathogenicity prediction 
- [MutScore](https://mutscore-wgt7hvakhq-ew.a.run.app/) - missense variant pathogenicity prediction 
- [SpliceAI](https://github.com/Illumina/SpliceAI) - deep learning method for predicting splice impacts
  - [SpliceAI-Lookup](https://spliceailookup.broadinstitute.org/) - easy SpliceAI lookup
- [ConSplice](https://github.com/mikecormier/ConSplice) - splice variant impact predictions modeled on constraint
- [DeepSEA](http://deepsea.princeton.edu/) - deep learning approach to predict variant chromatin and regulatory effects
- [COSMIS](https://github.com/CapraLab/cosmis) - 3D mutational constraint on amino acid sites in the human proteome
- [ExPecto](https://github.com/FunctionLab/ExPecto) - deep learning approach to predict variant expression effects in specific tissues
- [Basenji](https://github.com/calico/basenji) - deep CNN approach to predict variant chromatin and regulatory effects
- [VEST4](https://karchinlab.org/apps/appVest.html) - missense variant pathogenicity prediction

### literature and clinical knowledge 
- [LitVar](https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/LitVar/#!?query=) - variant based literature searching
- [GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1116/) - expert reviewed gene and genetic disease information
- [UpToDate](https://www.uptodate.com/contents/search) - clinical summaries and reviews of disease (subscription required)
- [Genomenon](https://www.genomenon.com/) - variant interpretation, associated literature, automatic notification if there are new publications
- [GHR](https://ghr.nlm.nih.gov/) - genetics home resource of genetic diseases written more for a lay audience

### specific data sets and visualizations
- [gnomAD](https://gnomad.broadinstitute.org/) - population scale variant allele frequency,loss of function intolerance scoring and structural variation
- [RGC](https://rgc-research.regeneron.com/) - million exome variant browser from regeneron
- [Pangolin](https://github.com/broadinstitute/pangolin) - deep-learning based method for predicting splice site strengths with scores for gnomAD
- [CCR](https://s3.us-east-2.amazonaws.com/ccrs/ccr.html) - constrained coding regions metric visual portal
- [DECIPHER](https://decipher.sanger.ac.uk/) - gene:disease association from Deciphering Developmental Disorders (DDD) Study
- [PERviewer](http://per.broadinstitute.org/) - pathogenic variant enriched regions across genes and gene families in gnomAD
- [BRAVO](https://bravo.sph.umich.edu/) - web tool for exploring TOPMed data
- [dbVar](https://www.ncbi.nlm.nih.gov/dbvar) - database of strutural variation
- [DGV](http://dgv.tcag.ca/dgv/app/home) - database of strutural variation
- [ExAC](http://exac.broadinstitute.org/) - superseeded by Gnomad
- [1000Genomes](http://www.internationalgenome.org/data/) - superseeded by Gnomad

### gene:phenotype information
- [HPO](https://hpo.jax.org/) - human phenotype ontology
- [GTR](https://www.ncbi.nlm.nih.gov/gtr/) - genetic testing registry, phenotype, disease, genetic test information
- [GWAS Catalog](https://www.ebi.ac.uk/gwas/) - catalog of GWAS associations from EBI/NHGRI
- [PheGenI](https://www.ncbi.nlm.nih.gov/gap/phegeni/) - gene:disease associations using underlying GWAS data
- [PhenoScanner](http://www.phenoscanner.medschl.cam.ac.uk/phenoscanner) - gene and variant disease/phenotype association using publicly available datasets
- [genebass](https://genebass.org/) - web browsing of UK Biobank genotype:phenotype data
- [AZPheWas](https://azphewas.com/) - gene:phenotype associations from the UK Biobank
- [EGA](https://ega-archive.org/) - aggregate genome:phenome associations
- [AMELIE](https://amelie.stanford.edu/) - Mendelian disease gene prioritization based on literature
- [Phenolyzer](http://phenolyzer.wglab.org/) - gene list generation using disease and phenotype terms
- [GeneNetwork](https://www.genenetwork.nl/) - gene networks from HPO terms
- [Mendelian](https://app.mendelian.co/) - phenotype based differential diagnosis tool including gene list generation
- [PubCaseFinder](https://pubcasefinder.dbcls.jp/) - phenotype based differential diagnosis and suspected disease predictor
- [PhenotypeGenerator](https://www.kimg.eu/generator/) - phenotype based gene list generation using HPO terms
- [PanelApp](https://panelapp.genomicsengland.co.uk/) - expert crowdsourced gene panels
- [gene2phenotype](https://www.ebi.ac.uk/gene2phenotype) - phenotype and biological process information using gene inputs
- [GADO](https://genenetwork.nl/gado/) - gene list generation and prioritization given a set of HPO terms
- [Phenomizer](http://compbio.charite.de/phenomizer/) - phenotype based gene list generation using HPO terms

### gene:drug interactions and pharmacogenomics
- [PAnno](https://github.com/PreMedKB/PAnno) - pharmacogenomic variant annotation and interpretation from a VCF input
- [PharmCat](https://github.com/PharmGKB/PharmCAT) - pharmacogenomic variant annotation and interpretation from a VCF input
- [Open Targets](https://genetics.opentargets.org/) - gene:disease/phenotype association focused on drug targets
- [PHAROS](https://pharos.ncats.nih.gov/) - gene:drug interaction and information, the druggable genome
- [TargetValidation](https://www.targetvalidation.org/) - gene and phenotype based genetic associations, pathways, drug targets
- [DGIdb](http://www.dgidb.org/) - drug:gene interaction database for potential therapeutic targets

### single cell data
- [CellxGene](https://cellxgene.cziscience.com/) - single cell data visualization and download across a large number of datasets
- [SingleCellPortal](https://singlecell.broadinstitute.org/single_cell) - data portal for finding single cell datasets across organisms and diseases
- [MouseBrain](http://mousebrain.org/) - single cell expression data through mouse prenatal and postnatal development
- [GutCellAtlas](https://www.gutcellatlas.org/) - single cell expression data across gut tissues and disease states
- [Descartes](https://descartes.brotmanbaty.org/) - single cell expression data across human development and chromatin state, plus other organisms
- [HumanFetalGutAtlas](https://simmonslab.shinyapps.io/FetalAtlasDataPortal/) - single cell expression data through human fetal gut development

### genetic and protein interactions
- [HumanBase](https://humanbase.io/) - machine learning driven biological knowledge
- [DIDA](http://dida.ibsquare.be/) - curated digenic disease database and gene pairs for oligogenic inheritance 
- [dbPTM](https://awi.cuhk.edu.cn/dbPTM/) - database of protein post-translational modifications

### gene and protein expression, pathway information
- [ProteinAtlas](http://www.proteinatlas.org/) - RNA and protein level gene expression across tissues and subcellular localization
- [GTEx](https://gtexportal.org/) - RNA and protein level gene expression across tissues
- [DICE](https://dice-database.org/) - gene expression and eQTLs in specific immune cell types
- [BioPlex](http://bioplex.hms.harvard.edu/) - protein:protein interaction networks database
- [InnateDB](http://www.innatedb.ca/) - protein:protein interaction networks database
- [BioGRID](https://thebiogrid.org/) - protein:protein interaction networks database
- [IntAct](https://www.ebi.ac.uk/intact/) - molecular interaction database
- [PathwayCommons](https://www.pathwaycommons.org/) - compiled pathway, interaction database
- [Reactome](https://reactome.org/) - signaling and regulatory pathway database with visualizations
- [KEGG](http://www.genome.jp/kegg/pathway.html) - biochemical pathway visualizations and information
- [signalink](http://signalink.org/) - signaling pathway cross-talks, transcription factors, miRNAs and regulatory enzymes
- [GO](http://geneontology.org/) - gene ontology, biological pathway enrichment, biological process, molecular function
- [Cytoscape](https://cytoscape.org/) - open source software for network and pathway visualization

### other

- [Normalizer](https://mutalyzer.nl/normalizer) - "normalize" a variant description to HGVS
- [PhenoTips](https://phenotips.org/) - for managing phenotype and pedigree information
- [OpenPedigree](https://github.com/phenotips/open-pedigree) - free open source pedigree drawing software
- [Visidata](https://www.visidata.org/) - visualize and interact with csv, tsv and other data in the terminal
- [TRUSST](https://www.grnpedia.org/trrust/) - transcriptional regulatory relationships using text mining
- [MyGene2](https://mygene2.org/MyGene2/) - gene and variant matching
- [GeneMatcher](https://genematcher.org/) - gene:disease matchmaking resource
- [ModelMatcher](https://www.modelmatcher.net/) - gene and model organism matching for variants/genes of uncertain significance
- [denovo-db](http://denovo-db.gs.washington.edu/denovo-db/) - database of de novo variants from publicly available datasets
- [VarMap](https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/DisaStr/GetPage.pl?varmap=TRUE) - genomic coordinate mapping to protein sequence and structure
- [RCSB-PDB](https://www.rcsb.org/) - protein structure database
- [BeaconNetwork](https://beacon-network.org/) - variant matching across numerous large databases, absent phenotype information
- [HUGO](https://genenames.org/) - accepted gene nomenclature and naming conventions
- [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) - liftover coordinates between genome builds
- [DataMed](https://datamed.org/) - searching across large published bioinformatics data sets
- [HSF3](http://www.umd.be/HSF3/index.html) - human splice analysis and predictor
- [Smart](https://smart.servier.com/) - free medical and biological art for presentations 
- [BioRender](https://biorender.com/) - biological figure and presentation creation

## command line tools

- [slivar](https://github.com/brentp/slivar) - rapid variant filtering and annotation for prioritization
- [VEP](https://uswest.ensembl.org/info/docs/tools/vep/script/index.html) - variant consequence and annotation, including plugins and custom plugins
- [vcfanno](https://github.com/brentp/vcfanno) - annotate vcfs with other vcfs/beds/bams
- [CADD-SV](https://github.com/kircherlab/CADD-SV) - pathogenicity scoring for SVs
- [SVAFotate](https://github.com/fakedrtom/SVAFotate) - structural variant population allele frequency annotation
- [manta](https://github.com/Illumina/manta) - sv calling
- [smoove](https://github.com/brentp/smoove) - structural variant calling, but smoothly
- [bcftools](https://github.com/samtools/bcftools) - vcf and bcf manipulation
- [samtools](https://github.com/samtools) - bam and cram manipulation
- [AnnotSV](https://lbgi.fr/AnnotSV/) - SV annotation into easy formats
- [ClassifyCNV](https://github.com/Genotek/ClassifyCNV) - classify CNVs according to ACMG criteria
- [RUFUS](https://github.com/jandrewrfarrell/RUFUS) - kmer based de novo variant caller
- [novoCaller](https://github.com/bgm-cwg/novoCaller) - Bayesian inspired de novo variant caller
- [PathoPredictor](https://github.com/samesense/pathopredictor) - missense variant classifier
- [HipSTR](https://github.com/tfwillems/HipSTR) - genotype short tandem repeats
- [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) - estimating repeat expansion sizes
- [GEMINI](https://github.com/arq5x/gemini) - builds a queryable database with comprehensive variant annotations
