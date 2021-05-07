//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.potential.groovy

import ffx.potential.cli.PotentialScript

import org.biojava.nbio.structure.StructureIO
import org.biojava.nbio.structure.Structure
import org.biojava.nbio.structure.Chain
import org.biojava.nbio.structure.Group
import org.biojava.nbio.structure.Atom
import org.biojava.nbio.structure.GroupType
import org.biojava.nbio.structure.PDBStatus

import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import java.util.concurrent.TimeUnit

/**
 * The ReadProteinNet reads a text-based protein net and outputs B-Factors
 *
 * <br>
 * Usage:
 * <br>
 * ffxc ReadProteinNet [options] &lt;filename&gt;
 */
@Command(description = " ReadProteinNet reads a text-based ProteinNet record and adds B-Factors to the dataset.", name = "ffxc ReadProteinNet")
class ReadProteinNet extends PotentialScript {

    /**
     * Argument should be one text-based ProteinNet file.
     */
    @Parameters(arity = "1", paramLabel = "file",
            description = 'Text-based ProteinNet file.')
    List<String> filenames = null

    /**
     * -d or --dummy 
     */
    @CommandLine.Option(names = ["-d", "--dummy"], paramLabel = "false",
            description = "Fill B-Factor arrays with 1.0 to act as dummy values for testing/validation evaluations")
    private boolean useDummy = false

    /**
     * --noNMR
     */
    @CommandLine.Option(names = ["--noNMR"], paramLabel = "false",
            description = "Leave NMR structures out of dataset")
    private boolean noNMR = false

    /**
     * --noNMR
     */
    @CommandLine.Option(names = ["--normalize"], paramLabel = "false",
            description = "Normalize B-Factors from 0 to 1 for relative uncertainty predictions.")
    private boolean normalize = false

    /**
    * --ids
    */
    @CommandLine.Option(names = ["--ids"], paramLabel = "file", description = "File with labels to include")
    private String labels = null


    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    String nextLine
    BufferedReader br

    int totalProteins = 0
    int badProteins = 0
    int notCrystals = 0
    int noCAs = 0
    int noPDBs = 0
    int noAlignments = 0
    int memoryErrors = 0
    int oneModels = 0
    List<String> noPDB
    List<String> noCA
    List<String> noAlignment
    List<String> notCrystal
    List<String> memoryError
    List<String> oneModel

    String currentID
    String currentSeq
    int seqLength
    int originalLength

    List<Integer> maskIdx

    HashMap<String, String> obsoleteMap
    HashSet<String> idSet

    /**
     * Execute the script.
     */
    @Override
    ReadProteinNet run() {

        if (!init()) {
            return this
        }

        String inputPath = filenames.get(0)

        logger.info("\n Opening ProteinNet " + inputPath)
        logger.info("")
        br = new BufferedReader(new FileReader(inputPath))

        noPDB = new ArrayList<>()
        noCA = new ArrayList<>()
        noAlignment = new ArrayList<>()
        notCrystal = new ArrayList<>()
        memoryError = new ArrayList<>()
        oneModel = new ArrayList<>()

        maskIdx = new ArrayList<>()

        String totalRecord = ""

        obsoleteMap = new HashMap<String, String>()

        File obsoleteFile = new File("/Dedicated/schnieders/gqi/folding/training/ProteinNetX/CASP12/training/obsolete.csv");
        Scanner myReader = new Scanner(obsoleteFile);
        while (myReader.hasNextLine()) {
            String data = myReader.nextLine();
            dataArr = data.split(",")
            if (dataArr.length > 1) {
                obsoleteMap.put(dataArr[0],dataArr[1])
            } else {
                obsoleteMap.put(dataArr[0],dataArr[0])
            }
        }

        if (labels != null) {
            idSet = new HashSet<String>()
            File labelFile = new File(labels)
            Scanner labelReader = new Scanner(labelFile)
            while (labelReader.hasNextLine()) {
                String data = labelReader.nextLine()
                idSet.add(data)
            }
        }
        loop:
        while (true) {

            String currentRecord = readRecord()

            if (useDummy) {
                if (currentRecord != null) {
                    totalProteins += 1
                    currentRecord = currentRecord.concat("[BFACTOR]\n")

                    double[] bFactorArrFinalN = new double[originalLength]
                    double[] bFactorArrFinalCA = new double[originalLength]
                    double[] bFactorArrFinalC = new double[originalLength]
                    Arrays.fill(bFactorArrFinalN, (double) 1.0)
                    Arrays.fill(bFactorArrFinalCA, (double) 1.0)
                    Arrays.fill(bFactorArrFinalC, (double) 1.0)

                    for (int i = 0; i < originalLength; i++) {
                        currentRecord = currentRecord.concat(Double.toString(bFactorArrFinalN[i]) + "\t")
                        currentRecord = currentRecord.concat(Double.toString(bFactorArrFinalCA[i]) + "\t")
                        currentRecord = currentRecord.concat(Double.toString(bFactorArrFinalC[i]) + "\t")
                    }

                    currentRecord = currentRecord.concat("\n\n")
                    totalRecord = totalRecord.concat(currentRecord)

                } else {
                    break loop
                }
            } else {

                if (currentRecord != null) {
                    totalProteins += 1
                    currentRecord = currentRecord.concat("[BFACTOR]\n")

                    List<Double> bFactorArrInitN = new ArrayList<>()
                    List<Double> bFactorArrInitCA = new ArrayList<>()
                    List<Double> bFactorArrInitC = new ArrayList<>()

                    String temp = ""

                    for (int i = 0; i < currentSeq.length(); i++) {
                        if (!maskIdx.contains(i)) {
                            temp += currentSeq.charAt(i)
                        }
                    }

                    currentSeq = temp
                    seqLength = currentSeq.length()

                    String bioJavaSeq = ""

                    try {
                        boolean write = true
                        String[] idArr = currentID.split("_")
                        String pdbID = idArr[0]
                        if (pdbID.contains("#")) {
                            pdbID = pdbID.split("#")[1]
                        }
                        Structure s1
                        currentPDBID = obsoleteMap.get(pdbID)
                        if (currentPDBID != null) {
                            s1 = StructureIO.getStructure(currentPDBID)
                        } else {
                            s1 = StructureIO.getStructure(pdbID)
                        }

                        int numResidues = 0

                        int numModels = s1.nrModels()

                        if (s1.isCrystallographic()) {
                            for (Chain chain : s1.getChains()) {
                                for (Group group : chain.getAtomGroups(GroupType.AMINOACID)) {
                                    Atom n = group.getAtom("N")
                                    Atom ca = group.getAtom("CA")
                                    Atom c = group.getAtom("C")
                                    if (n != null && ca != null && c != null) {
                                        bFactorArrInitN.add(n.getTempFactor())
                                        bFactorArrInitCA.add(ca.getTempFactor())
                                        bFactorArrInitC.add(c.getTempFactor())
                                    } else {
                                        bFactorArrInitN.add(0.0)
                                        bFactorArrInitCA.add(0.0)
                                        bFactorArrInitC.add(0.0)
                                    }
                                    bioJavaSeq += group.getAminoType()
                                    numResidues += 1
                                }
                            }
                        } else if (s1.isNmr()) {
                            if (!noNMR) {
                                List<Double[]> nList = new ArrayList<>()
                                List<Double[]> caList = new ArrayList<>()
                                List<Double[]> cList = new ArrayList<>()

                                if (numModels > 1) {
                                    for (int i = 0; i < numModels; i++) {
                                        for (Chain chain : s1.getModel(i)) {
                                            int atomCounter = 0
                                            for (Group group : chain.getAtomGroups(GroupType.AMINOACID)) {
                                                Atom n = group.getAtom("N")
                                                Atom ca = group.getAtom("CA")
                                                Atom c = group.getAtom("C")
                                                if (i == 0) {
                                                    if (n != null && ca != null && c != null) {
                                                        nList.add(n.getCoords())
                                                        caList.add(ca.getCoords())
                                                        cList.add(c.getCoords())
                                                    } else {
                                                        nList.add([0.0, 0.0, 0.0])
                                                        caList.add([0.0, 0.0, 0.0])
                                                        cList.add([0.0, 0.0, 0.0])
                                                    }
                                                    bioJavaSeq += group.getAminoType()
                                                    numResidues += 1
                                                } else {
                                                    if (n != null && ca != null && c != null) {
                                                        for (int j = 0; j < 3; j++) {
                                                            nList.get(atomCounter)[j] += n.getCoords()[j]
                                                            caList.get(atomCounter)[j] += ca.getCoords()[j]
                                                            cList.get(atomCounter)[j] += c.getCoords()[j]
                                                        }
                                                    }
                                                }
                                                atomCounter++
                                            }
                                        }
                                    }

                                    for (int i = 0; i < numResidues; i++) {
                                        for (int j = 0; j < 3; j++) {
                                            nList.get(i)[j] = nList.get(i)[j] / numModels
                                            caList.get(i)[j] = caList.get(i)[j] / numModels
                                            cList.get(i)[j] = cList.get(i)[j] / numModels
                                        }
                                    }

                                    List<Double> nSumList = new ArrayList<>()
                                    List<Double> caSumList = new ArrayList<>()
                                    List<Double> cSumList = new ArrayList<>()

                                    for (int i = 0; i < numModels; i++) {
                                        for (Chain chain : s1.getModel(i)) {
                                            int atomCounter = 0
                                            for (Group group : chain.getAtomGroups(GroupType.AMINOACID)) {
                                                Atom n = group.getAtom("N")
                                                Atom ca = group.getAtom("CA")
                                                Atom c = group.getAtom("C")
                                                double[] diff_sq = new double[3]
                                                for (int j = 0; j < 3; j++) {
                                                    if (n != null) {
                                                        diff_sq[j] = Math.pow(n.getCoords()[j] - nList.get(atomCounter)[j], 2)
                                                    } else {
                                                        diff_sq[j] = 0.0
                                                    }
                                                }
                                                if (i == 0) {
                                                    nSumList.add(diff_sq[0] + diff_sq[1] + diff_sq[2])
                                                } else {
                                                    nSumList.set(atomCounter, nSumList.get(atomCounter) + diff_sq[0] + diff_sq[1] + diff_sq[2])
                                                }

                                                diff_sq = new double[3]
                                                for (int j = 0; j < 3; j++) {
                                                    if (ca != null) {
                                                        diff_sq[j] = Math.pow(ca.getCoords()[j] - caList.get(atomCounter)[j], 2)
                                                    } else {
                                                        diff_sq[j] = 0.0
                                                    }
                                                }
                                                if (i == 0) {
                                                    caSumList.add(diff_sq[0] + diff_sq[1] + diff_sq[2])
                                                } else {
                                                    caSumList.set(atomCounter, caSumList.get(atomCounter) + diff_sq[0] + diff_sq[1] + diff_sq[2])
                                                }

                                                diff_sq = new double[3]
                                                for (int j = 0; j < 3; j++) {
                                                    if (c != null) {
                                                        diff_sq[j] = Math.pow(c.getCoords()[j] - cList.get(atomCounter)[j], 2)
                                                    } else {
                                                        diff_sq[j] = 0.0
                                                    }
                                                }
                                                if (i == 0) {
                                                    cSumList.add(diff_sq[0] + diff_sq[1] + diff_sq[2])
                                                } else {
                                                    cSumList.set(atomCounter, cSumList.get(atomCounter) + diff_sq[0] + diff_sq[1] + diff_sq[2])
                                                }
                                                atomCounter++
                                            }
                                        }
                                    }

                                    double bFactorConstant = (8.0 * Math.pow(Math.PI, 2)) / 3.0

                                    for (int i = 0; i < numResidues; i++) {
                                        bFactorArrInitN.add(bFactorConstant * (nSumList.get(i) / numModels) + 20.0)
                                        bFactorArrInitCA.add(bFactorConstant * (caSumList.get(i) / numModels) + 20.0)
                                        bFactorArrInitC.add(bFactorConstant * (cSumList.get(i) / numModels) + 20.0)
                                    }

                                } else {
                                    oneModel.add(currentID)
                                    oneModels += 1
                                    badProteins += 1
                                    write = false

                                }
                            } else {
                                notCrystal.add(currentID)
                                notCrystals += 1                                
                                badProteins += 1
                                write = false
                            }      
                        } else {
                            notCrystal.add(currentID)
                            notCrystals += 1
                            badProteins += 1
                            write = false                            
                        }

                        double[] bFactorArrFinalN = new double[seqLength]
                        double[] bFactorArrFinalCA = new double[seqLength]
                        double[] bFactorArrFinalC = new double[seqLength]
                        Arrays.fill(bFactorArrFinalN, (double) -1.0)
                        Arrays.fill(bFactorArrFinalCA, (double) -1.0)
                        Arrays.fill(bFactorArrFinalC, (double) -1.0)

                        if (write && currentSeq.equals(bioJavaSeq)) {

                            for (int i = 0; i < seqLength; i++) {
                                bFactorArrFinalN[i] = (double) bFactorArrInitN.get(i)
                                bFactorArrFinalCA[i] = (double) bFactorArrInitCA.get(i)
                                bFactorArrFinalC[i] = (double) bFactorArrInitC.get(i)

                                if (bFactorArrFinalN[i] < 1e-10 ||
                                        bFactorArrFinalCA[i] < 1e-10 ||
                                        bFactorArrFinalC[i] < 1e-10) {
                                    badProteins += 1
                                    noCA.add(currentID)
                                    noCAs += 1
                                    write = false
                                    break
                                }

                                if (Math.abs(bFactorArrFinalN[i]) < 1.0) {
                                    bFactorArrFinalN[i] = (double) 1.0
                                }

                                if (Math.abs(bFactorArrFinalCA[i]) < 1.0) {
                                    bFactorArrFinalCA[i] = (double) 1.0
                                }

                                if (Math.abs(bFactorArrFinalC[i]) < 1.0) {
                                    bFactorArrFinalC[i] = (double) 1.0
                                }

                            }
                        } else if (write && !(currentSeq.equals(bioJavaSeq))) {
                            if (bioJavaSeq.length() == 0) {
                                badProteins += 1
                                noAlignments += 1
                                noAlignment.add(currentID + " " + seqLength.toString() + " " + numResidues.toString())
                                write = false
                            } else {
                                if (seqLength > numResidues) {
                                    int startIndex = currentSeq.indexOf(bioJavaSeq)
                                    if (startIndex != -1) {
                                        for (int i = 0; i < numResidues; i++) {
                                            bFactorArrFinalN[i + startIndex] = (double) bFactorArrInitN.get(i)
                                            bFactorArrFinalCA[i + startIndex] = (double) bFactorArrInitCA.get(i)
                                            bFactorArrFinalC[i + startIndex] = (double) bFactorArrInitC.get(i)

                                            if (bFactorArrFinalN[i + startIndex] < 1e-10 ||
                                                    bFactorArrFinalCA[i + startIndex] < 1e-10 ||
                                                    bFactorArrFinalC[i + startIndex] < 1e-10) {
                                                badProteins += 1
                                                noCA.add(currentID)
                                                noCAs += 1
                                                write = false
                                                break
                                            }
                                            if (Math.abs(bFactorArrFinalN[i + startIndex]) < 1.0) {
                                                bFactorArrFinalN[i + startIndex] = (double) 1.0
                                            }

                                            if (Math.abs(bFactorArrFinalCA[i + startIndex]) < 1.0) {
                                                bFactorArrFinalCA[i + startIndex] = (double) 1.0
                                            }

                                            if (Math.abs(bFactorArrFinalC[i + startIndex]) < 1.0) {
                                                bFactorArrFinalC[i + startIndex] = (double) 1.0
                                            }                                       
                                        }
                                    } else {
                                        badProteins += 1
                                        noAlignments += 1
                                        noAlignment.add(currentID + " " + seqLength.toString() + " " + numResidues.toString())
                                        write = false
                                    }
                                } else {
                                    int startIndex = bioJavaSeq.indexOf(currentSeq)
                                    if (startIndex != -1) {
                                        for (int i = 0; i < seqLength; i++) {
                                            bFactorArrFinalN[i] = (double) bFactorArrInitN.get(i + startIndex)
                                            bFactorArrFinalCA[i] = (double) bFactorArrInitCA.get(i + startIndex)
                                            bFactorArrFinalC[i] = (double) bFactorArrInitC.get(i + startIndex)

                                            if (bFactorArrFinalN[i] < 1e-10 ||
                                                    bFactorArrFinalCA[i] < 1e-10 ||
                                                    bFactorArrFinalC[i] < 1e-10) {
                                                badProteins += 1
                                                noCA.add(currentID)
                                                noCAs += 1
                                                write = false
                                                break
                                            }

                                            if (Math.abs(bFactorArrFinalN[i]) < 1.0) {
                                                bFactorArrFinalN[i] = (double) 1.0
                                            }

                                            if (Math.abs(bFactorArrFinalCA[i]) < 1.0) {
                                                bFactorArrFinalCA[i] = (double) 1.0
                                            }

                                            if (Math.abs(bFactorArrFinalC[i]) < 1.0) {
                                                bFactorArrFinalC[i] = (double) 1.0
                                            }

                                        }

                                    } else {
                                        badProteins += 1
                                        noAlignments += 1
                                        noAlignment.add(currentID + " " + seqLength.toString() + " " + numResidues.toString())
                                        write = false
                                    }
                                }
                            }
                        }

                        if (normalize) {
                            // Find the mean
                            double totalSum = 0
                            int count = 0
                            for (int i = 0; i < bFactorArrFinalN.length; i++) {
                                totalSum += bFactorArrFinalN[i]
                                count++
                            }
                            for (int i = 0; i < bFactorArrFinalCA.length; i++) {
                                totalSum += bFactorArrFinalCA[i]
                                count++
                            }
                            for (int i = 0; i < bFactorArrFinalC.length; i++) {
                                totalSum += bFactorArrFinalC[i]
                                count++
                            }
                            double mean = totalSum / count

                            // Find the standard deviation
                            double totalSquaredError = 0
                            for (int i = 0; i < bFactorArrFinalN.length; i++) {
                                double diff = bFactorArrFinalN[i] - mean
                                totalSquaredError += diff * diff
                            }
                            for (int i = 0; i < bFactorArrFinalCA.length; i++) {
                                double diff = bFactorArrFinalN[i] - mean
                                totalSquaredError += diff * diff
                            }
                            for (int i = 0; i < bFactorArrFinalC.length; i++) {
                                double diff = bFactorArrFinalN[i] - mean
                                totalSquaredError += diff * diff
                            }
                            double std = Math.sqrt(totalSquaredError / (count - 1))

                            // Normalize values using min/max rescaling
                            for (int i = 0; i < bFactorArrFinalN.length; i++) {
                                bFactorArrFinalN[i] = (bFactorArrFinalN[i] - mean) / std;
                            }
                            for (int i = 0; i < bFactorArrFinalCA.length; i++) {
                                bFactorArrFinalCA[i] = (bFactorArrFinalCA[i] - mean) / std;
                            }
                            for (int i = 0; i < bFactorArrFinalC.length; i++) {
                                bFactorArrFinalC[i] = (bFactorArrFinalC[i] - mean) / std;
                            }
                        }

                        if (write) {
                            int index = 0
                            for (int i = 0; i < originalLength; i++) {
                                if (maskIdx.contains(i)) {
                                    currentRecord = currentRecord.concat("-1.0\t-1.0\t-1.0\t")
                                } else {
                                    currentRecord = currentRecord.concat(Double.toString(bFactorArrFinalN[index]) + "\t")
                                    currentRecord = currentRecord.concat(Double.toString(bFactorArrFinalCA[index]) + "\t")
                                    currentRecord = currentRecord.concat(Double.toString(bFactorArrFinalC[index]) + "\t")
                                    index++
                                }
                            }
                            if (index != seqLength) {
                                println("we have a problem")
                                println(currentID)
                            }
                            currentRecord = currentRecord.concat("\n\n")
                            totalRecord = totalRecord.concat(currentRecord)
                        }                       

                    } catch (FileNotFoundException e) {
                        println("PDB not found for: " + currentID)
                        badProteins += 1
                        noPDB.add(currentID)
                        noPDBs += 1
                    } catch (OutOfMemoryError e) {
                        println("Out of memory: " + currentID)
                        badProteins += 1
                    } catch (Exception e) {
                        e.printStackTrace()
                        println(currentID)
                        badProteins += 1
                        noPDB.add(currentID)
                        noPDBs += 1
                    }
                } else {
                    break loop
                }
            }

        }

        File outfile = new File(inputPath + ".xnet")
        FileWriter fileWriter = new FileWriter(outfile)

        fileWriter.write(totalRecord)
        fileWriter.close()

        logger.info(" Total structures in ProteinNet: " + totalProteins)
        logger.info(" Total errors: " + badProteins)
        double percent = badProteins / totalProteins * 100
        logger.info(" Percent Error: " + percent)
        logger.info("")

        logger.info(" No PDB found corresponding to PBD ID: " + noPDBs)
        logger.info(" " + noPDB.toString())
        logger.info("")

        logger.info(" No alignment between primary sequence given and PDB sequence: " + noAlignments)
        logger.info(" " + noAlignment.toString())
        logger.info("")

        logger.info(" PDB does not correspond to a crystallographic or NMR structure: " + notCrystals)
        logger.info(" " + notCrystal.toString())
        logger.info("")

        logger.info(" NMR structure only has one model: " + oneModels)
        logger.info(" " + oneModel.toString())
        logger.info("")

        logger.info(" Java memory errors: " + memoryErrors)
        logger.info(" " + memoryError.toString())
        logger.info("")

        logger.info(" No C-alpha: " + noCAs)
        logger.info(" " + noCA.toString())
        logger.info("")

        return this
    }

    /**
     *
     * @return: String containing record for the protein
     */
    private String readRecord() {
        String record = ""
        while ((nextLine = br.readLine()) != null) {
            switch (nextLine) {
                case "[ID]":
                    record = record.concat(nextLine + "\n")
                    nextLine = br.readLine()
                    currentID = nextLine
                    record = record.concat(nextLine + "\n")
                    break;
                case "[PRIMARY]":
                    record = record.concat(nextLine + "\n")
                    nextLine = br.readLine()
                    currentSeq = nextLine
                    originalLength = nextLine.length()
                    record = record.concat(nextLine + "\n")
                    break;
                case "[MASK]":
                    maskIdx.clear()
                    record = record.concat(nextLine + "\n")
                    nextLine = br.readLine()
                    currentMask = nextLine
                    for (int i = 0; i < currentMask.length(); i++) {
                        Character c = currentMask.charAt(i)
                        if (c == '-') {
                            maskIdx.add(i)
                        }
                    }
                    record = record.concat(nextLine + "\n")
                    break;
                case "":
                    if (idSet != null) {
                        if (idSet.contains(currentID)) {
                            return record
                        } else {
                            record = ""
                            break;
                        }  
                    } else {
                        return record
                    }           
                default:
                    record = record.concat(nextLine + "\n")
                    break;
            }
        }
        return null
    }
}
