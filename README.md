# variant_prioritization
purpose: a curated list of resources and tools for variant prioritization in clinical and research settings, welcoming additions and changes

## web tools

- [MARRVEL](http://marrvel.org/) - model organism focused gene information
- [Varsome](https://varsome.com/) - aggregate gene and variant information 
- [HumanMine](http://www.humanmine.org/) - aggregate gene information
- [GeneCards](http://www.genecards.org/) - aggregate gene information
- [ClinGen](https://www.clinicalgenome.org/) - expert curated gene:disease associations
- [GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1116/) - expert reviewed gene and genetic disease information
- [MonarchInitiative](https://monarchinitiative.org/) - gene, disease, biological function information, including gene ortholog phenotypes
- [UpToDate](https://www.uptodate.com/contents/search) - clinical summaries and reviews of disease (subscription required)
- [Gnomad](https://gnomad.broadinstitute.org/) - population scale variant allele frequency,loss of function intolerance scoring and structural variation
- [OMIM](http://omim.org/) - online Mendelian inheritance in man, gene:disease association, gene function, variant molecular genetics
- [MedGen](https://www.ncbi.nlm.nih.gov/medgen/) - clinically focused gene and disease information including literature review
- [GenicIntolerance](http://genic-intolerance.org/) - gene based loss of function intolerance scoring
- [Phenomizer](http://compbio.charite.de/phenomizer/) - phenotype based gene list generation using HPO terms
- [Phenolyzer](http://phenolyzer.wglab.org/) - gene list generation using disease and phenotype terms
- [Mendelian](https://app.mendelian.co/) - phenotype based differential diagnosis tool including gene list generation
- [PubCaseFinder](https://pubcasefinder.dbcls.jp/) - phenotype based differential diagnosis and suspected disease predictor
- [PhenotypeGenerator](https://www.kimg.eu/generator/) - phenotype based gene list generation using HPO terms
- [HPO](https://hpo.jax.org/) - human phenotype ontology
- [GTR](https://www.ncbi.nlm.nih.gov/gtr/) - genetic testing registry, phenotype, disease, genetic test information
- [PanelApp](https://panelapp.genomicsengland.co.uk/) - expert crowdsourced gene panels
- [gene2phenotype](https://www.ebi.ac.uk/gene2phenotype) - phenotype and biological process information using gene inputs
- [PheGenI](https://www.ncbi.nlm.nih.gov/gap/phegeni/) - gene:disease associations using underlying GWAS data
- [PhenoScanner](http://www.phenoscanner.medschl.cam.ac.uk/phenoscanner) - gene and variant disease/phenotype association using publicly available datasets
- [Genomenon](https://www.genomenon.com/) - variant interpretation, associated literature, automatic notification if there are new publications
- [CCR](https://s3.us-east-2.amazonaws.com/ccrs/ccr.html) - constrained coding regions metric visual portal
- [ProteinAtlas](http://www.proteinatlas.org/) - RNA and protein level gene expression across tissues and subcellular localization
- [GTEx](https://gtexportal.org/) - RNA and protein level gene expression across tissues
- [BioPlex](http://bioplex.hms.harvard.edu/) - protein:protein interaction networks database
- [InnateDB](http://www.innatedb.ca/) - protein:protein interaction networks database
- [BioGRID](https://thebiogrid.org/) - protein:protein interaction networks database
- [TRUSST](https://www.grnpedia.org/trrust/) - transcriptional regulatory relationships using text mining
- [TargetValidation](https://www.targetvalidation.org/) - gene and phenotype based genetic associations, pathways, drug targets
- [DGIdb](http://www.dgidb.org/) - drug:gene interaction database for potential therapeutic targets
- [SNP](https://www.ncbi.nlm.nih.gov/snp/) - variant curation and information, dbsnp, rsID matching
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) - clinical variant interpretation across multiple clinical sites and commercial testing providers
- [DECIPHER](https://decipher.sanger.ac.uk/) - gene:disease association from Deciphering Developmental Disorders (DDD) Study
- [KEGG](http://www.genome.jp/kegg/pathway.html) - biochemical pathway visualizations and information
- [GO](http://geneontology.org/) - gene ontology, biological pathway enrichment, biological process, molecular function
- [signalink](http://signalink.org/) - signaling pathway cross-talks, transcription factors, miRNAs and regulatory enzymes
- [Reactome](https://reactome.org/) - signaling and regulatory pathway database with visualizations
- [denovo-db](http://denovo-db.gs.washington.edu/denovo-db/) - database of de novo variants from publicly available datasets
- [RCSB-PDB](https://www.rcsb.org/) - protein structure database
- [GeneMatcher](https://genematcher.org/) - gene:disease matchmaking resource
- [dbVar](https://www.ncbi.nlm.nih.gov/dbvar) - database of strutural variation
- [DGV](http://dgv.tcag.ca/dgv/app/home) - database of strutural variation
- [EBI](http://www.ebi.ac.uk/) - aggregate gene information, similar to NIH resources
- [BeaconNetwork](https://beacon-network.org/) - variant matching across numerous large databases, absent phenotype information
- [PhenIX](http://compbio.charite.de/PhenIX/) - vcf based variant prioritization considering allele frequency and phenotype information
- [GAVIN](https://molgenis20.gcc.rug.nl/menu/main/home) - gene aware variant interpretation using vcf input
- [HUGO](https://genenames.org/) - accepted gene nomenclature and naming conventions
- [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) - liftover coordinates between genome builds
- [ExAC](http://exac.broadinstitute.org/) - likely superseeded by Gnomad
- [1000Genomes](http://www.internationalgenome.org/data/) - likely superseeded by Gnomad
- [UCSCbrowser](https://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=724544049_pngh3ffiA9LYDPiWojHaNAcDu3CA) - supremely customizable genome browser
- [DataMed](https://datamed.org/) - searching across large published bioinformatics data sets
- [WashUStLEpigenome](http://epigenomegateway.wustl.edu/browser/) - epigenome browser and visualization
- [HSF3](http://www.umd.be/HSF3/index.html) - human splice analysis and predictor
- [GHR](https://ghr.nlm.nih.gov/) - genetics home resource of genetic diseases written more for a lay audience

## command line tools
#### including a typical, but not authoritative, usage

- [GEMINI](https://github.com/arq5x/gemini) - builds a queryable database with comprehensive variant annotations
```
gemini load -t VEP -v $VCF.vep.vcf.gz --cores 24 -p $PED $VCF.gem.db
```
- [slivar](https://github.com/brentp/slivar) - rapid variant filtering and annotation for prioritization
```
slivar_static expr --vcf $VCF --ped $PED --pass-only --js $JS -g $GNOMAD --trio "denovo:denovo(kid, mom, dad) && hqrv(variant, INFO, '0.01')" --trio "x_denovo:x_denovo(kid, mom, dad) && hqrv(variant, INFO, '0.01') && (variant.CHROM == 'X' || variant.CHROM == 'chrX')" --trio "recessive:recessive(kid, mom, dad) && hqrv(variant, INFO, '0.01')" --trio "x_recessive:x_recessive(kid, mom, dad) && hqrv(variant, INFO, '0.01') && (variant.CHROM == 'X' || variant.CHROM == 'chrX')"
```
