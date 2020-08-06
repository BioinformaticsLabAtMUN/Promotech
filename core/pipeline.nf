params.bacteria    = 'ECOLI_2';
params.outDir      = "$baseDir/outDir/${params.bacteria}"
params.indexTrain  = "$baseDir/train.py"

//Variables
indexBedFile       = file("$baseDir/bacteria/${params.bacteria}/${params.bacteria}.bed")
indexFastaFile     = file("$baseDir/bacteria/${params.bacteria}/${params.bacteria}.fasta")
lengthFile         = file("${params.bacteria}.genome")
trainFile          = file(params.indexTrain)

Channel
   .from(indexFastaFile)
   .splitFasta( record: [id: true, seqString: true ])
   .collectFile(name: lengthFile) { record -> record.id + "\t" + record.seqString.length() + "\n"}
   .set{indexGenomeFile1}

Channel
   .from(indexFastaFile)
   .splitFasta( record: [id: true, seqString: true ])
   .collectFile(name: lengthFile) { record -> record.id + "\t" + record.seqString.length() + "\n"}
   .set{indexGenomeFile2}


println """\
PROMOTER PREDICTION = ${params.bacteria}
===============================

AUTHOR: RUBEN CHEVEZ
DATE:   OCT 15, 2018

"""
.stripIndent()

process getPositivePromoterBED {
    // container 'genomicpariscentre/bedtools'
    input:
        file indexBedFile
        file indexFastaFile     
        file "${params.bacteria}.genome" from indexGenomeFile1
        publishDir params.outDir, mode: 'copy' 
    output:
        file "positive.bed" into (promoter_bed_file, promoter_bed_file2 )
    script: 
    """
    slopBed -i ${indexBedFile} -s -g ${params.bacteria}.genome -l 39 -r 0 > "positive.bed"
    """
}

process getPositivePromoterFASTA {
    // container 'genomicpariscentre/bedtools'
    input:
        file promoter_bed_file
        file indexFastaFile   
        publishDir params.outDir, mode: 'copy' 
    output:
        file "positive.fasta" into positive_promoter_sequence
    script:
    """
    bedtools getfasta -s -bed ${promoter_bed_file} -fi ${indexFastaFile}  -fo  "positive.fasta" 
    """
}

process getNegativePromotersBED {
    // container 'genomicpariscentre/bedtools'
    input:
        file promoter_bed_file2
        file positive_promoter_sequence 
        file "${params.bacteria}.genome" from indexGenomeFile2
        publishDir params.outDir, mode: 'copy' 
    output:
        file "negative.bed" into negative_promoter_bed
    script: 
    """
    bedtools random -l 40 -n `grep -i -v ">" ${positive_promoter_sequence} | wc -l` -g ${params.bacteria}.genome > "tmp_negative_promoter.bed"
    bedtools subtract -f 0.13 -s -A -a tmp_negative_promoter.bed -b ${promoter_bed_file2} > negative.bed
    """
}

process getNegativePromoterFASTA {
    // container 'genomicpariscentre/bedtools'
    input:
        file negative_promoter_bed
        file indexFastaFile   
        publishDir params.outDir, mode: 'copy' 
    output:
        file "negative.fasta" into negative_promoter_sequence
    script:
    """
    bedtools getfasta -s  -bed ${negative_promoter_bed} -fi ${indexFastaFile}  -fo  "negative.fasta" 
    """
}