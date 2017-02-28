package edu.caltech.lncrna.clap.programs;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import edu.caltech.lncrna.bio.alignment.PairedEndAlignment;
import edu.caltech.lncrna.bio.annotation.BedFileRecord;
import edu.caltech.lncrna.bio.annotation.Block;
import edu.caltech.lncrna.bio.datastructures.GenomeTree;
import edu.caltech.lncrna.bio.io.BedParser;
import edu.caltech.lncrna.bio.io.PairedEndBamParser;

public class EnrichmentDistribution {

    private Path inputPath;
    private Path exptPath;
    private Path genesPath;
    private Path posOutputPath;
    private Path negOutputPath;
    private int tileLength;
    private int stepSize;
    private long inputNorm;
    private long exptNorm;
    private GenomeTree<PairedEndAlignment> inputReads;
    private GenomeTree<PairedEndAlignment> exptReads;
    
    private static final int MINIMUM_NUMBER_READS = 10;

    public static void main(String[] args) throws IOException {
        (new EnrichmentDistribution())
            .parseArgs(args)
            .loadReads()
            .calculateValuesAndPrint();
	}

    private EnrichmentDistribution parseArgs(String[] args) {
        Options options = new Options();

        options.addOption(Option.builder()
                .longOpt("input")
                .desc("input BAM file")
                .argName("BAM")
                .hasArg()
                .required()
                .build());

        options.addOption(Option.builder()
                .longOpt("expt")
                .desc("experimental BAM file")
                .argName("BAM")
                .hasArg()
                .required()
                .build());

        options.addOption(Option.builder()
                .longOpt("genes")
                .desc("BED file of genes or transcripts")
                .argName("BED")
                .hasArg()
                .required()
                .build());

        options.addOption(Option.builder()
                .longOpt("tile")
                .desc("length of each tile")
                .argName("INT")
                .hasArg()
                .required()
                .build());

        options.addOption(Option.builder()
                .longOpt("step")
                .desc("step size for tiling")
                .argName("INT")
                .hasArg()
                .required()
                .build());

        options.addOption(Option.builder()
                .longOpt("inputnorm")
                .desc("input normalization factor")
                .argName("INT")
                .hasArg()
                .required()
                .build());
        
        options.addOption(Option.builder()
                .longOpt("exptnorm")
                .desc("expt normalization factor")
                .argName("INT")
                .hasArg()
                .required()
                .build());
        
        options.addOption(Option.builder()
                .longOpt("pos")
                .desc("positive output file")
                .argName("OUT")
                .hasArg()
                .required()
                .build());
        
        options.addOption(Option.builder()
                .longOpt("neg")
                .desc("negative output file")
                .argName("OUT")
                .hasArg()
                .required()
                .build());

        try {
            CommandLineParser parser = new DefaultParser();
            CommandLine cmd = parser.parse(options, args);
            inputPath = Paths.get(cmd.getOptionValue("input"));
            exptPath = Paths.get(cmd.getOptionValue("expt"));
            genesPath = Paths.get(cmd.getOptionValue("genes"));
            posOutputPath = Paths.get(cmd.getOptionValue("pos"));
            negOutputPath = Paths.get(cmd.getOptionValue("neg"));
            tileLength = Integer.parseInt(cmd.getOptionValue("tile"));
            stepSize = Integer.parseInt(cmd.getOptionValue("step"));
            inputNorm = Long.parseLong(cmd.getOptionValue("inputnorm"));
            exptNorm = Long.parseLong(cmd.getOptionValue("exptnorm"));
        } catch (ParseException | NumberFormatException e) {
            printHelp(options);
            System.exit(-1);
        }
        
        return this;
    }
    
    private void printHelp(Options opts) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.setDescPadding(0);
        String header = "\n";
        String footer = "\n";
        formatter.printHelp("java -jar EnrichmentDistribution.jar", header,
                opts, footer, true);
    }

    /**
     * Loads reads from input and elution paired-end BAM files.
     * 
     * Parses the reads from the input and elution BAM files, filters out the
     * unmapped read-pairs, and store the remaining mapped read-pairs in
     * memory. 
     */
    private EnrichmentDistribution loadReads() {
        inputReads = new GenomeTree<PairedEndAlignment>();
        try (PairedEndBamParser p = new PairedEndBamParser(inputPath)) {
            p.stream()
             .map(x -> x.getAlignment())
             .filter(x -> x.isPresent())
             .map(x -> x.get())
             .forEach(x -> inputReads.insert(x));
        }

        exptReads = new GenomeTree<PairedEndAlignment>();
        try (PairedEndBamParser p = new PairedEndBamParser(exptPath)) {
            p.stream()
             .map(x -> x.getAlignment())
             .filter(x -> x.isPresent())
             .map(x -> x.get())
             .forEach(x -> exptReads.insert(x));
        }
        
        return this;
    }

    /**
     * Calculates enrichment values and prints to output.
     * 
     * This method creates a tiling for each gene in the input BED file.
     * For each tile, this method calculates read-enrichment over that tile;
     * then, depending on the strandedness of the gene, prints the enrichment
     * value to either a positive or negative output file.
     *
     * @throws IOException
     */
    private EnrichmentDistribution calculateValuesAndPrint()
            throws IOException {

        try (BedParser p = new BedParser(genesPath);
             BufferedWriter posWriter =
                     Files.newBufferedWriter(posOutputPath);
             BufferedWriter negWriter =
                     Files.newBufferedWriter(negOutputPath)) {

            while (p.hasNext()) {
                BedFileRecord b = p.next();
                calculateValueOverBedRecordAndPrint(b, posWriter, negWriter);
            }
        }
        return this;
    }
    
    /**
     * Calculates and prints enrichment values over a gene represented in a
     * single BED file record.
     * 
     * This calculation is done for the tilings of each exon in the gene.
     * @throws IOException
     */
    private void calculateValueOverBedRecordAndPrint(BedFileRecord bed,
            BufferedWriter posOut, BufferedWriter negOut) throws IOException {

        Iterator<Block> exons = bed.getBlockIterator();

        while (exons.hasNext()) {
            Block exon = exons.next();
            calculateValueOverExonAndPrint(exon, posOut, negOut);
        }
    }
    
    /**
     * Tiles an exon, calculates enrichment values over that tiling, and
     * prints to output.
     * @throws IOException
     */
    private void calculateValueOverExonAndPrint(Block exon,
            BufferedWriter posOut, BufferedWriter negOut) throws IOException {
        
        Iterator<Block> tiles = exon.tile(tileLength, stepSize);
        
        while (tiles.hasNext()) {
            Block tile = tiles.next();
            calculateValueOverTileAndPrint(tile, posOut, negOut);
        }
    }
    
    /**
     * Calculates the enrichment value over a single tile and prints to output.
     */
    private void calculateValueOverTileAndPrint(Block tile,
            BufferedWriter posOut, BufferedWriter negOut) throws IOException {
 
        switch (tile.getStrand()) {
        case POSITIVE:
            calculateValueAndPrint(tile, posOut);
            break;
        case NEGATIVE:
            calculateValueAndPrint(tile, negOut);
            break;
        default:
            // do nothing
        }
    }
    
    /**
     * Calculates and prints enrichment value to the given writer.
     */
    private void calculateValueAndPrint(Block tile, BufferedWriter writer)
            throws IOException {
        
        long numInputReads = inputReads.getNumOverlappers(tile);

        if (numInputReads >= MINIMUM_NUMBER_READS) {
            long numExptReads = exptReads.getNumOverlappers(tile);
            int midpoint = (tile.getStart() + tile.getEnd()) / 2;
            writer.write(tile.getReferenceName());
            writer.write("\t");
            writer.write(String.valueOf(midpoint));
            writer.write("\t");
            writer.write(String.valueOf(midpoint + 1));
            writer.write("\t");
            writer.write(String.valueOf(((double) numExptReads * inputNorm) /
                    (numInputReads * exptNorm)));
            writer.write("\t");
            writer.write(String.valueOf((double) numInputReads / inputNorm));
            writer.write("\t");
            writer.write(String.valueOf(numInputReads));
            writer.write("\t");
            writer.write(String.valueOf(numExptReads));
            writer.newLine();
        }
    }
}