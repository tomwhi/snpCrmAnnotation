-- MySQL dump 10.13  Distrib 5.5.44, for debian-linux-gnu (x86_64)
--
-- Host: localhost    Database: motifs
-- ------------------------------------------------------
-- Server version	5.5.44-0ubuntu0.14.04.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `DBD_BlockItems`
--

DROP TABLE IF EXISTS `DBD_BlockItems`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `DBD_BlockItems` (
  `pfamID` varchar(255) DEFAULT NULL,
  `blockID` int(11) DEFAULT NULL,
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `start` int(11) DEFAULT NULL,
  `end` int(11) DEFAULT NULL,
  `eValue` double DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `fk_DBD_BlockItems_1` (`pfamID`),
  KEY `fk_DBD_BlockItems_2` (`blockID`),
  CONSTRAINT `fk_DBD_BlockItems_1` FOREIGN KEY (`pfamID`) REFERENCES `PFAM_DBDs` (`pfamID`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `fk_DBD_BlockItems_2` FOREIGN KEY (`blockID`) REFERENCES `dbdBlocks` (`id`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=66962 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `PFAM_DBDs`
--

DROP TABLE IF EXISTS `PFAM_DBDs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `PFAM_DBDs` (
  `pfamID` varchar(255) NOT NULL,
  `Name` varchar(255) DEFAULT NULL,
  `Description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`pfamID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aaSeq`
--

DROP TABLE IF EXISTS `aaSeq`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `aaSeq` (
  `seqID` int(11) NOT NULL AUTO_INCREMENT,
  `aaSeq` varchar(767) DEFAULT NULL,
  PRIMARY KEY (`seqID`),
  UNIQUE KEY `aaSeq_UNIQUE` (`aaSeq`)
) ENGINE=InnoDB AUTO_INCREMENT=16067 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aaSeqSim`
--

DROP TABLE IF EXISTS `aaSeqSim`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `aaSeqSim` (
  `seq1id` int(11) NOT NULL,
  `seq2id` int(11) NOT NULL,
  `identity` double DEFAULT NULL,
  `similarity` double DEFAULT NULL,
  PRIMARY KEY (`seq1id`,`seq2id`),
  KEY `fk_dbdSim_1` (`seq1id`),
  KEY `fk_dbdSim_2` (`seq2id`),
  CONSTRAINT `fk_dbdSim_1` FOREIGN KEY (`seq1id`) REFERENCES `aaSeq` (`seqID`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `fk_dbdSim_2` FOREIGN KEY (`seq2id`) REFERENCES `aaSeq` (`seqID`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `affyProbes`
--

DROP TABLE IF EXISTS `affyProbes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `affyProbes` (
  `GeneID` varchar(255) NOT NULL,
  `probeset` varchar(255) NOT NULL,
  PRIMARY KEY (`GeneID`,`probeset`),
  KEY `probeset` (`probeset`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `breastEnrichedProbesets`
--

DROP TABLE IF EXISTS `breastEnrichedProbesets`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `breastEnrichedProbesets` (
  `probeset` varchar(255) NOT NULL DEFAULT '',
  PRIMARY KEY (`probeset`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `dbdBlocks`
--

DROP TABLE IF EXISTS `dbdBlocks`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `dbdBlocks` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `proteinID` varchar(255) DEFAULT NULL,
  `aaSeqID` int(11) DEFAULT NULL,
  `start` int(11) DEFAULT NULL,
  `end` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `fk_dbdSeq_1` (`proteinID`),
  KEY `fk_dbdBlocks_1` (`aaSeqID`),
  CONSTRAINT `fk_dbdBlocks_1` FOREIGN KEY (`aaSeqID`) REFERENCES `aaSeq` (`seqID`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `fk_dbdSeq_1` FOREIGN KEY (`proteinID`) REFERENCES `proteins` (`proteinID`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=32176 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `geneProteins`
--

DROP TABLE IF EXISTS `geneProteins`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `geneProteins` (
  `geneID` varchar(255) NOT NULL,
  `proteinID` varchar(255) NOT NULL,
  PRIMARY KEY (`geneID`,`proteinID`),
  KEY `fk_geneProteins_1` (`geneID`),
  KEY `fk_geneProteins_2` (`proteinID`),
  CONSTRAINT `fk_geneProteins_1` FOREIGN KEY (`geneID`) REFERENCES `genes` (`GeneID`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `fk_geneProteins_2` FOREIGN KEY (`proteinID`) REFERENCES `proteins` (`proteinID`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `genes`
--

DROP TABLE IF EXISTS `genes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genes` (
  `GeneID` varchar(255) NOT NULL,
  `Name` varchar(255) DEFAULT NULL,
  `Species` varchar(45) DEFAULT NULL,
  PRIMARY KEY (`GeneID`),
  KEY `nameIdx` (`Name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gnf`
--

DROP TABLE IF EXISTS `gnf`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gnf` (
  `probeset` varchar(255) NOT NULL,
  `tissue` varchar(255) NOT NULL,
  `expnLevel` double DEFAULT NULL,
  PRIMARY KEY (`probeset`,`tissue`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `kgAlias`
--

DROP TABLE IF EXISTS `kgAlias`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `kgAlias` (
  `kgID` varchar(40) NOT NULL DEFAULT '',
  `alias` varchar(80) DEFAULT NULL,
  KEY `kgID` (`kgID`),
  KEY `alias` (`alias`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `kgXref`
--

DROP TABLE IF EXISTS `kgXref`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `kgXref` (
  `kgID` varchar(255) NOT NULL,
  `mRNA` varchar(255) NOT NULL,
  `spID` varchar(255) NOT NULL,
  `spDisplayID` varchar(255) NOT NULL,
  `geneSymbol` varchar(255) NOT NULL,
  `refseq` varchar(255) NOT NULL,
  `protAcc` varchar(255) NOT NULL,
  `description` longblob NOT NULL,
  `rfamAcc` varchar(255) NOT NULL,
  `tRnaName` varchar(255) NOT NULL,
  KEY `kgID` (`kgID`),
  KEY `mRNA` (`mRNA`),
  KEY `spID` (`spID`),
  KEY `spDisplayID` (`spDisplayID`),
  KEY `geneSymbol` (`geneSymbol`),
  KEY `refseq` (`refseq`),
  KEY `protAcc` (`protAcc`),
  KEY `rfamAcc` (`rfamAcc`),
  KEY `tRnaName` (`tRnaName`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `knownGene`
--

DROP TABLE IF EXISTS `knownGene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `knownGene` (
  `name` varchar(255) NOT NULL DEFAULT '',
  `chrom` varchar(255) NOT NULL DEFAULT '',
  `strand` char(1) NOT NULL DEFAULT '',
  `txStart` int(10) unsigned NOT NULL DEFAULT '0',
  `txEnd` int(10) unsigned NOT NULL DEFAULT '0',
  `cdsStart` int(10) unsigned NOT NULL DEFAULT '0',
  `cdsEnd` int(10) unsigned NOT NULL DEFAULT '0',
  `exonCount` int(10) unsigned NOT NULL DEFAULT '0',
  `exonStarts` longblob NOT NULL,
  `exonEnds` longblob NOT NULL,
  `proteinID` varchar(40) NOT NULL DEFAULT '',
  `alignID` varchar(255) NOT NULL DEFAULT '',
  KEY `name` (`name`),
  KEY `chrom` (`chrom`(16),`txStart`),
  KEY `chrom_2` (`chrom`(16),`txEnd`),
  KEY `protein` (`proteinID`(16)),
  KEY `align` (`alignID`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `motifParams`
--

DROP TABLE IF EXISTS `motifParams`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `motifParams` (
  `rc` tinyint(1) NOT NULL,
  `trimIC` double NOT NULL,
  `motifPrefix` varchar(255) NOT NULL,
  PRIMARY KEY (`rc`,`trimIC`,`motifPrefix`),
  KEY `fk_motifParams_1` (`motifPrefix`),
  CONSTRAINT `fk_motifParams_1` FOREIGN KEY (`motifPrefix`) REFERENCES `motifs` (`prefix`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `motifs`
--

DROP TABLE IF EXISTS `motifs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `motifs` (
  `prefix` varchar(255) NOT NULL,
  `type` varchar(255) DEFAULT NULL,
  `heterodimer` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`prefix`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `probeset2gene`
--

DROP TABLE IF EXISTS `probeset2gene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `probeset2gene` (
  `probeset` varchar(255) NOT NULL,
  `geneName` varchar(255) NOT NULL,
  PRIMARY KEY (`probeset`,`geneName`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `proteins`
--

DROP TABLE IF EXISTS `proteins`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `proteins` (
  `proteinID` varchar(255) NOT NULL,
  `aaSeq` varchar(65000) DEFAULT NULL,
  PRIMARY KEY (`proteinID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `teemusReducedSELEX`
--

DROP TABLE IF EXISTS `teemusReducedSELEX`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `teemusReducedSELEX` (
  `motifPrefix` varchar(255) NOT NULL,
  PRIMARY KEY (`motifPrefix`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tfMotifs`
--

DROP TABLE IF EXISTS `tfMotifs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tfMotifs` (
  `motifPrefix` varchar(255) NOT NULL,
  `geneID` varchar(255) NOT NULL,
  PRIMARY KEY (`motifPrefix`,`geneID`),
  KEY `fk_tfMotifs_2` (`geneID`),
  KEY `fk_tfMotifs_1` (`motifPrefix`),
  CONSTRAINT `fk_tfMotifs_1` FOREIGN KEY (`motifPrefix`) REFERENCES `motifs` (`prefix`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `fk_tfMotifs_2` FOREIGN KEY (`geneID`) REFERENCES `genes` (`GeneID`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmp1`
--

DROP TABLE IF EXISTS `tmp1`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tmp1` (
  `Name` varchar(255) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmp2`
--

DROP TABLE IF EXISTS `tmp2`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tmp2` (
  `Name` varchar(255) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmp2_jaspar_identity100`
--

DROP TABLE IF EXISTS `tmp2_jaspar_identity100`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tmp2_jaspar_identity100` (
  `geneName` varchar(255) DEFAULT NULL,
  `domainName` varchar(255) DEFAULT NULL,
  `prefix` varchar(255) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmp2_jaspar_sim90`
--

DROP TABLE IF EXISTS `tmp2_jaspar_sim90`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tmp2_jaspar_sim90` (
  `geneName` varchar(255) DEFAULT NULL,
  `domainName` varchar(255) DEFAULT NULL,
  `prefix` varchar(255) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmp2_selex_identity100`
--

DROP TABLE IF EXISTS `tmp2_selex_identity100`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tmp2_selex_identity100` (
  `geneName` varchar(255) DEFAULT NULL,
  `domainName` varchar(255) DEFAULT NULL,
  `prefix` varchar(255) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmp2_selex_sim90`
--

DROP TABLE IF EXISTS `tmp2_selex_sim90`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tmp2_selex_sim90` (
  `geneName` varchar(255) DEFAULT NULL,
  `domainName` varchar(255) DEFAULT NULL,
  `prefix` varchar(255) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmp2_uniprobe_identity100`
--

DROP TABLE IF EXISTS `tmp2_uniprobe_identity100`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tmp2_uniprobe_identity100` (
  `geneName` varchar(255) DEFAULT NULL,
  `domainName` varchar(255) DEFAULT NULL,
  `prefix` varchar(255) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmp2_uniprobe_sim90`
--

DROP TABLE IF EXISTS `tmp2_uniprobe_sim90`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tmp2_uniprobe_sim90` (
  `geneName` varchar(255) DEFAULT NULL,
  `domainName` varchar(255) DEFAULT NULL,
  `prefix` varchar(255) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmpGeneIDs`
--

DROP TABLE IF EXISTS `tmpGeneIDs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tmpGeneIDs` (
  `geneID` varchar(255) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmpGeneIDs2`
--

DROP TABLE IF EXISTS `tmpGeneIDs2`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tmpGeneIDs2` (
  `GeneID` varchar(255) NOT NULL,
  `similarity` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmpGenesWithDBD`
--

DROP TABLE IF EXISTS `tmpGenesWithDBD`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tmpGenesWithDBD` (
  `GeneID` varchar(255) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmpTFs`
--

DROP TABLE IF EXISTS `tmpTFs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tmpTFs` (
  `GeneID` varchar(255) NOT NULL,
  `Name` varchar(255) DEFAULT NULL,
  `Species` varchar(45) DEFAULT NULL,
  `similarity` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `u133a_tmp`
--

DROP TABLE IF EXISTS `u133a_tmp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `u133a_tmp` (
  `id` varchar(255) NOT NULL DEFAULT '',
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `u133a_tmp2`
--

DROP TABLE IF EXISTS `u133a_tmp2`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `u133a_tmp2` (
  `name` varchar(255) NOT NULL DEFAULT '',
  PRIMARY KEY (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `vaqTF_DBDs`
--

DROP TABLE IF EXISTS `vaqTF_DBDs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `vaqTF_DBDs` (
  `EnsembleID` varchar(255) NOT NULL,
  `InterproID` varchar(255) NOT NULL,
  PRIMARY KEY (`EnsembleID`,`InterproID`),
  KEY `fk_vaqTF_DBDs_1` (`EnsembleID`),
  KEY `fk_vaqTF_DBDs_2` (`InterproID`),
  CONSTRAINT `fk_vaqTF_DBDs_1` FOREIGN KEY (`EnsembleID`) REFERENCES `vaqTFs` (`EnsembleID`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `fk_vaqTF_DBDs_2` FOREIGN KEY (`InterproID`) REFERENCES `vaq_DBD_IP` (`InterproID`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `vaqTFs`
--

DROP TABLE IF EXISTS `vaqTFs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `vaqTFs` (
  `EnsembleID` varchar(255) NOT NULL,
  `confidence` varchar(45) DEFAULT NULL,
  PRIMARY KEY (`EnsembleID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `vaq_DBD_IP`
--

DROP TABLE IF EXISTS `vaq_DBD_IP`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `vaq_DBD_IP` (
  `InterproID` varchar(255) NOT NULL,
  `Type` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`InterproID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `vaq_IPI`
--

DROP TABLE IF EXISTS `vaq_IPI`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `vaq_IPI` (
  `IPI_ID` varchar(255) NOT NULL,
  `EnsembleID` varchar(255) NOT NULL,
  PRIMARY KEY (`IPI_ID`,`EnsembleID`),
  KEY `fk_vaq_IPI_1` (`EnsembleID`),
  CONSTRAINT `fk_vaq_IPI_1` FOREIGN KEY (`EnsembleID`) REFERENCES `vaqTFs` (`EnsembleID`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `vaq_InterproAnnot`
--

DROP TABLE IF EXISTS `vaq_InterproAnnot`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `vaq_InterproAnnot` (
  `InterproID` varchar(255) NOT NULL,
  `Family` varchar(255) DEFAULT NULL,
  `Description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`InterproID`),
  KEY `fk_vaq_InterproAnnot_1` (`InterproID`),
  CONSTRAINT `fk_vaq_InterproAnnot_1` FOREIGN KEY (`InterproID`) REFERENCES `vaq_DBD_IP` (`InterproID`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2016-04-14 15:29:24
