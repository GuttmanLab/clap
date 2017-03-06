package edu.caltech.lncrna.clap.programs;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import edu.caltech.lncrna.bio.alignment.PairedEndAlignment;
import edu.caltech.lncrna.bio.annotation.Block;
import edu.caltech.lncrna.bio.annotation.Populated;
import edu.caltech.lncrna.bio.annotation.Strand;
import edu.caltech.lncrna.bio.annotation.WindowIterator;
import edu.caltech.lncrna.bio.io.PairedEndBamParser;

public class WindowedEnrichmentCalculator {

    private Path inputPath;
    private Path elutionPath;
    private Path posOutputPath;
    private Path negOutputPath;
    private int tileLength;
    private int stepSize;
    private final Map<Block, Integer> inputWindows = new HashMap<>();
    
    public static void main(String[] args) throws IOException {
        (new WindowedEnrichmentCalculator())
            .parseArgs(args)
            .calculateEnrichments(Strand.POSITIVE)
            .calculateEnrichments(Strand.NEGATIVE);
	}

    private WindowedEnrichmentCalculator parseArgs(String[] args) {
        Options options = new Options();

        options.addOption(Option.builder()
                .longOpt("input")
                .desc("input BAM file")
                .argName("BAM")
                .hasArg()
                .required()
                .build());

        options.addOption(Option.builder()
                .longOpt("elution")
                .desc("elution BAM file")
                .argName("BAM")
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
            elutionPath = Paths.get(cmd.getOptionValue("elution"));
            posOutputPath = Paths.get(cmd.getOptionValue("pos"));
            negOutputPath = Paths.get(cmd.getOptionValue("neg"));
            tileLength = Integer.parseInt(cmd.getOptionValue("tile"));
            stepSize = Integer.parseInt(cmd.getOptionValue("step"));
            
            if (tileLength < 1 || stepSize < 1) {
                throw new NumberFormatException();
            }
            
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
        formatter.printHelp("java -jar WindowedEnrichmentCalculator.jar", header,
                opts, footer, true);
    }

    private WindowedEnrichmentCalculator calculateEnrichments(Strand strand)
            throws IOException {

        calculateInputCoverage(strand);
        calculateElutionCoverageAndPrintEnrichment(strand);
        return this;
    }
    
    private void calculateInputCoverage(Strand strand) {
        try (PairedEndBamParser readPairs = 
                new PairedEndBamParser(inputPath)) {
   
            Iterator<PairedEndAlignment> alignments =
                    readPairs.getAlignmentStream()
                             .filter(x -> x.getStrand().equals(strand))
                             .iterator();
            
            WindowIterator<PairedEndAlignment> windows =
                    new WindowIterator<>(alignments, tileLength, stepSize);
            
            while (windows.hasNext()) {
                Populated<PairedEndAlignment> window = windows.next();
                inputWindows.put(new Block(window), window.getPopulationSize());
            }
        }
    }

    private void calculateElutionCoverageAndPrintEnrichment(Strand strand)
            throws IOException {

        try (BufferedWriter writer = strand.equals(Strand.POSITIVE)
                 ? Files.newBufferedWriter(posOutputPath)
                 : Files.newBufferedWriter(negOutputPath)) {
            
            try (PairedEndBamParser readPairs = new PairedEndBamParser(elutionPath)) {
            
                Iterator<PairedEndAlignment> alignments =
                        readPairs.getAlignmentStream()
                                 .filter(x -> x.getStrand().equals(strand))
                                 .iterator();
            
                WindowIterator<PairedEndAlignment> windows =
                        new WindowIterator<>(alignments, tileLength, stepSize);
            
                while (windows.hasNext()) {
                    calculateElutionCoverageForOneWindowAndPrint(windows.next(),
                            writer);
                }
            }
            
            printRemainingInputWindows(writer);
        }
    }
    
    private void calculateElutionCoverageForOneWindowAndPrint(
            Populated<PairedEndAlignment> window, BufferedWriter writer)
            throws IOException {

        Block block = new Block(window);
        Integer matchingInputScore = inputWindows.remove(block);
        if (matchingInputScore == null) {
            matchingInputScore = 0;
        }
        
        writer.write(window.getReferenceName() + "\t" +
                     window.getStart() + "\t" +
                     window.getEnd() + "\t" +
                     matchingInputScore + "\t" +
                     window.getPopulationSize());
        writer.newLine();
    }
    
    private void printRemainingInputWindows(BufferedWriter writer)
            throws IOException {

        for (Map.Entry<Block, Integer> entry : inputWindows.entrySet()) {
            Block window = entry.getKey();
            int inputScore = entry.getValue();
            writer.write(window.getReferenceName() + "\t" +
                         window.getStart() + "\t" +
                         window.getEnd() + "\t" +
                         inputScore + "\t" +
                         0);
            writer.newLine();
        }
        
        inputWindows.clear();
    }
}