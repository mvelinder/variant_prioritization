# variant_prioritization
purpose: a curated list (in no particular order) of resources and tools for variant prioritization in clinical and research settings, always welcoming additions and changes

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
  - [CNV Interpretation](http://cnvcalc.clinicalgenome.org/cnvcalc/) - CNV pathogenicity calculator
- [HGMD](http://www.hgmd.cf.ac.uk/ac/index.php) - curated gene:disease associations based on literature, with public and paid version
- [MedGen](https://www.ncbi.nlm.nih.gov/medgen/) - clinically focused gene and disease information including literature review
- [SimpleClinVar](http://simple-clinvar.broadinstitute.org/) - simplied summaries and views of ClinVar data

### biological impact predictions
- [EVE](https://evemodel.org/) - evolutionary model of missense variant effects
- [REVEL](https://sites.google.com/site/revelgenomics/) - missense variant pathogenicity prediction 
- [COSMIS](https://github.com/CapraLab/cosmis) - 3D mutational constraint on amino acid sites in the human proteome
- [MutScore](https://mutscore-wgt7hvakhq-ew.a.run.app/) - missense variant pathogenicity prediction 
- [SpliceAI](https://github.com/Illumina/SpliceAI) - deep learning method for predicting splice impacts
  - [SpliceAI Lookup](https://spliceailookup.broadinstitute.org/) - easy SpliceAI lookup
- [ConSplice](https://github.com/mikecormier/ConSplice) - splice variant impact predictions modeled on constraint
- [DeepSEA](http://deepsea.princeton.edu/) - deep learning approach to predict variant chromatin and regulatory effects
- [ExPecto](https://github.com/FunctionLab/ExPecto) - deep learning approach to predict variant expression effects in specific tissues
- [Basenji](https://github.com/calico/basenji) - deep CNN approach to predict variant chromatin and regulatory effects
- [VEST4](https://karchinlab.org/apps/appVest.html) - missense variant pathogenicity prediction

### literature and clinical knowledge 
- [LitVar](https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/LitVar/#!?query=) - variant based literature searching
- [GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1116/) - expert reviewed gene and genetic disease information
- [UpToDate](https://www.uptodate.com/contents/search) - clinical summaries and reviews of disease (subscription required)
- [Genomenon](https://www.genomenon.com/) - variant interpretation, associated literature, automatic notification if there are new publications
- [GHR](https://ghr.nlm.nih.gov/) - genetics home resource of genetic diseases written more for a lay audience

### specific data sets and visualuzations
- [gnomAD](https://gnomad.broadinstitute.org/) - population scale variant allele frequency,loss of function intolerance scoring and structural variation
- [CCR](https://s3.us-east-2.amazonaws.com/ccrs/ccr.html) - constrained coding regions metric visual portal
- [DECIPHER](https://decipher.sanger.ac.uk/) - gene:disease association from Deciphering Developmental Disorders (DDD) Study
- [genebass](https://genebass.org/) - web browsing of UK Biobank genotype:phenotype data
- [AZPheWas](https://azphewas.com/) - gene:phenotype associations from the UK Biobank
- [EGA](https://ega-archive.org/) - aggregate genome:phenome associations
- [PERviewer](http://per.broadinstitute.org/) - pathogenic variant enriched regions across genes and gene families in gnomAD
- [PheGenI](https://www.ncbi.nlm.nih.gov/gap/phegeni/) - gene:disease associations using underlying GWAS data
- [PhenoScanner](http://www.phenoscanner.medschl.cam.ac.uk/phenoscanner) - gene and variant disease/phenotype association using publicly available datasets
- [BRAVO](https://bravo.sph.umich.edu/) - web tool for exploring TOPMed data
- [dbVar](https://www.ncbi.nlm.nih.gov/dbvar) - database of strutural variation
- [DGV](http://dgv.tcag.ca/dgv/app/home) - database of strutural variation
- [ExAC](http://exac.broadinstitute.org/) - superseeded by Gnomad
- [1000Genomes](http://www.internationalgenome.org/data/) - superseeded by Gnomad

### gene:phenotype information
- [HPO](https://hpo.jax.org/) - human phenotype ontology
- [GTR](https://www.ncbi.nlm.nih.gov/gtr/) - genetic testing registry, phenotype, disease, genetic test information
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

### gene:drug interactions
- [PHAROS](https://pharos.ncats.nih.gov/) - gene:drug interaction and information, the druggable genome
- [TargetValidation](https://www.targetvalidation.org/) - gene and phenotype based genetic associations, pathways, drug targets
- [DGIdb](http://www.dgidb.org/) - drug:gene interaction database for potential therapeutic targets

### genetic and protein interactions
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

### biological effect computational 

- [PhenoTips](https://phenotips.org/) - free open source software for managing phenotype and pedigree information 
- [TRUSST](https://www.grnpedia.org/trrust/) - transcriptional regulatory relationships using text mining
- [MyGene2](https://mygene2.org/MyGene2/) - gene and variant matching
- [GeneMatcher](https://genematcher.org/) - gene:disease matchmaking resource
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
