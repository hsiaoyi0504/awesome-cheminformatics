# Awesome Cheminformatics [![Awesome](https://awesome.re/badge.svg)](https://awesome.re)

> Cheminformatics (also known as chemoinformatics, chemioinformatics and chemical informatics) is the use of computer and informational techniques applied to a range of problems in the field of chemistry.— [Wikipedia](https://en.wikipedia.org/wiki/Cheminformatics)

A curated list of awesome Cheminformatics software, resources, and libraries. Mostly command line based, and free or open-source. Please feel free to [contribute](CONTRIBUTING.md) !

__tags:__
* Open Source: :book:
* Free: :credit_card:
* Comercial: :moneybag:

## Contents

* [Applications](#applications)
  * [Visualization](#visualization)
  * [Quantum Chemistry](#quantum)
  * [Toolkits](#toolkits)
  * [Standarize](#format)
  * [Molecular Descriptors](#des)
  * [Machine Learning](#ml)
  * [Databases](#db)
  * [Docking](#docking)
  * [Molecular Dynamics](#md)
  * [Protein](#protein)
* [Journals](#journals)
* [Resources](#resources)
  * [Courses](#courses)
  * [Blogs](#blogs)
  * [Books](#books)
* [See Also](#see-also)

## Applications

<a id="visualization"></a>

### Visualization

#### 3D-Visulization
* [PyMOL](https://sourceforge.net/projects/pymol/) - Python-enhanced molecular graphics tool. :book: :moneybag:
* [Jmol](http://jmol.sourceforge.net/) - Browser-based HTML5 viewer and stand-alone Java viewer for chemical structures in 3D. :book:
* [VMD](http://www.ks.uiuc.edu/Research/vmd/) - Molecular visualization program for displaying, animating, and analyzing large biomolecular systems using 3-D graphics and built-in scripting. :book:
* [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/) - Highly extensible program for interactive molecular visualization and analysis. :book:
* [3Dmol.js](https://3dmol.csb.pitt.edu/) - A modern, object-oriented JavaScript library for visualizing molecular data. :book:

#### 2D-Visulization
* [JSME](https://www.peter-ertl.com/jsme/) - JSME is a free molecule editor written in JavaScript. JSME is a direct successor of the JME Molecule Editor applet. :book:
* [ChemDoddle](https://www.chemdoodle.com/) - ChemDoddle is a cross-platform molecule editor. :moneybag:
* [MarvinJS](https://chemaxon.com/products/marvin-js) - Marvin JS provides quick and convenient ways to draw and modify standard and advanced chemical structures. It's seamlessly integrated into third-party web-based applications, and runs smoothly on all major browsers. It is a ChemAxon product. :credit_card: :moneybag:
* [Kekule.js](http://partridgejiang.github.io/Kekule.js/) - Front-end JavaScript library for providing the ability to represent, draw, edit, compare and search molecule structures on web browsers. :book:
* [JChemPaint](https://github.com/JChemPaint/jchempaint) - Chemical 2D structure editor application/applet based on the [Chemistry Development Kit](https://sourceforge.net/projects/cdk/). :book:
* [rdeditor](https://github.com/EBjerrum/rdeditor) - Simple RDKit molecule editor GUI using PySide. :book:

<a id="quantum"></a>

### Quantum Chemistry

* [Gaussian](http://gaussian.com/) - The most famous commericial quantum chemistry computational engine :moneybag:
* [GAMMES(US)](https://www.msg.chem.iastate.edu/gamess/) - The General Atomic and Molecular Electronic Structure System (GAMESS) is a general ab initio quantum chemistry package. :book:
* [Juagar](https://www.schrodinger.com/jaguar) - Rapid ab initio electronic structure package. :moneybag:
* [Q-Chem](http://q-chem.com/) - Q-Chem is a comprehensive ab initio quantum chemistry package for accurate predictions of molecular structures, reactivities, and vibrational, electronic and NMR spectra. :moneybag:
* [PSI4](http://www.psicode.org/) - PSI4 is an open-source suite of ab initio quantum chemistry programs designed for efficient, high-accuracy simulations of a variety of molecular properties. It is very easy to use and has an optional Python interface. :book:
* [cclib](https://github.com/cclib/cclib) - Parsers and algorithms for computational chemistry logfiles. 

<a id="toolkits"></a>

### Toolkits

* [RDKit](http://www.rdkit.org/) - Collection of cheminformatics and machine-learning software written in C++ and Python. :book:
* [Indigo](https://github.com/epam/Indigo) - Universal molecular toolkit that can be used for molecular fingerprinting, substructure search, and molecular visualization written in C++ package, with Java, C#, and Python wrappers. :book:
* [CDK (Chemistry Development Kit)](https://cdk.github.io) - Algorithms for structural chemo- and bioinformatics, implemented in Java. :book:
* [ChemmineR](https://www.bioconductor.org/packages/release/bioc/vignettes/ChemmineR/inst/doc/ChemmineR.html) - Cheminformatics package for analyzing drug-like small molecule data in R. :book:
* [Open Babel](http://openbabel.org/wiki/Main_Page) - Chemical toolbox designed to speak the many languages of chemical data. :book:
* [MayaChemTools](http://www.mayachemtools.org/index.html) - Collection of Perl and Python scripts, modules, and classes that support day-to-day computational discovery needs. :book:

<a id="format"></a>

### Standarize

* [standardiser](https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/) - Tool designed to provide a simple way of standardising molecules as a prelude to e.g. molecular modelling exercises. Github [address](https://github.com/flatkinson/standardiser). :book:
* [MolVS](https://github.com/mcs07/MolVS) - Molecule validation and standardization based on [RDKit](http://www.rdkit.org/). :book:
* [rd_filters](https://github.com/PatWalters/rd_filters) - This script provides a simple means of applying the functional group filters from the ChEMBL database, as well as a number of property filters from the RDKit, to a set of compounds. :book:


<a id="des"></a>

### Molecular Descriptors

* [mordred](https://github.com/mordred-descriptor/mordred) - Molecular descriptor calculator based on [RDKit](http://www.rdkit.org/). :book:
* [mol2vec](https://github.com/samoturk/mol2vec) - Vector representations of molecular substructures. :book:
* [Align-it](http://silicos-it.be.s3-website-eu-west-1.amazonaws.com/software/align-it/1.0.4/align-it.html#alignit-generating-pharmacophore-points) - Align molecules according their pharmacophores. :book:
* [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) - A software to calculate molecular descriptors and fingerprints. The software currently calculates 1875 descriptors (1444 1D, 2D descriptors and 431 3D descriptors) and 12 types of fingerprints (total 16092 bits). :book:

<a id="ml"></a>

### Machine Learning

* [DeepChem](https://github.com/deepchem/deepchem) - Democratizing deep learning for Chemistry. :book:
* [Tensor-Mol](https://github.com/jparkhill/TensorMol) - Tensor-Mol is a neural net work model to do moleuclar dynamics and energy estimation. :book:
* [keras-moelcules](https://github.com/maxhodak/keras-molecules) - A Keras implementation of Aspuru-Guzik's molecular autoencoder paper. :book:
* [molecular-autoencoder](https://github.com/HIPS/molecule-autoencoder) - Automatic chemical design using a data-driven continuous representation of molecules. [arxiv.org](https://arxiv.org/abs/1610.02415). :book:
* [SMILES-enumeration](https://github.com/EBjerrum/SMILES-enumeration) - SMILES enumeration is the process of writing out all possible SMILES forms of a molecule. It's a useful technique for data augmentation before sequence based modeling of molecules. [arxiv.org](https://arxiv.org/abs/1703.07076). :book:
* [grammar-VAE](https://github.com/mkusner/grammarVAE) - Code for the "Grammar Variational Autoencoder" [arxiv.org](https://arxiv.org/abs/1703.01925). :book:
* [ChemTS](https://github.com/tsudalab/ChemTS) - Molecule Design using Monte Carlo Tree Search with Neural Rollout. ChemTS can design novel molecules with desired properties(such as, HOMO-LUMO gap, energy, logp..). Combining with rDock, ChemTS can design molecules active to target proteins. [arxiv.org](https://arxiv.org/abs/1710.00616). :book:
* [reLeaSE](https://github.com/isayev/ReLeaSE) - Deep Reinforcement Learning for de-novo Drug Design. [arxiv.org](http://dx.doi.org/10.1126/sciadv.aap7885). :book:
* [chemprop](https://github.com/swansonk14/chemprop) - This repository contains message passing neural networks for molecular property prediction. :book:

<a id="db"></a>

### Databases

#### Databases

* [PubChem](https://pubchem.ncbi.nlm.nih.gov/) - PubChem is world's largest collection of freely accessible chemical information. Search chemicals by name, molecular formula, structure, and other identifiers. Find chemical and physical properties, biological activities, safety and toxicity information, patents, literature citations and more. :book:
* [ChEMBL](https://www.ebi.ac.uk/chembl/) - ChEMBL is a manually curated database of bioactive molecules with drug-like properties. It brings together chemical, bioactivity and genomic data to aid the translation of genomic information into effective new drugs. :book:
* [Excape](https://sandbox.ideaconsult.net/search/excape) - ExCAPE-DB is a large public chemogenomics dataset based on the PubChem and ChEMBL databases, and large scale standardisation (including tautomerization) of chemical structures using open source cheminformatics software. :book:
* [ZINC](http://zinc.docking.org/) - ZINC is a free database of commercially-available compounds for virtual screening.   ZINC contains over 35 million purchasable compounds in ready-to-dock, 3D formats.   ZINC is provided by the Irwin and Shoichet Laboratories in the Department of Pharmaceutical Chemistry at the University of California, San Francisco (UCSF). :book:
* [ZINC15](http://zinc15.docking.org/) - Welcome to ZINC, a free database of commercially-available compounds for virtual screening. ZINC contains over 230 million purchasable compounds in ready-to-dock, 3D formats. ZINC also contains over 750 million purchasable compounds you can search for analogs in under a minute.  :book:
* [ChemSpider](http://www.chemspider.com/) - ChemSpider is a free chemical structure database providing fast access to over 34 million structures, properties and associated information. :book:
* [DrugBank](https://www.drugbank.ca/) - The DrugBank database is a unique bioinformatics and cheminformatics resource that combines detailed drug data with comprehensive drug target information. Acadamic free. :moneybag:
* [KEGG](https://www.genome.jp/kegg/) - KEGG is a database resource for understanding high-level functions and utilities of the biological system, such as the cell, the organism and the ecosystem, from molecular-level information, especially large-scale molecular datasets generated by genome sequencing and other high-throughput experimental technologies. It also provide [drug database](https://www.genome.jp/kegg/drug/) and [compound database](https://www.genome.jp/kegg/compound/). :book:

#### Web APIs

* [webchem](https://github.com/ropensci/webchem) - Chemical Information from the Web. :book:
* [PubChemPy](http://pubchempy.readthedocs.io) - Python wrapper for the PubChem PUG REST API. :book:
* [ChemSpiPy](http://chemspipy.readthedocs.org) - Python wrapper for the ChemSpider API. :book:
* [CIRpy](http://cirpy.readthedocs.org/) - Python wrapper for the NCI Chemical Identifier Resolver (CIR). :book:
* [Beaker](https://github.com/chembl/chembl_beaker) - [RDKit](http://www.rdkit.org/) and [OSRA](https://cactus.nci.nih.gov/osra/) in the [Bottle](http://bottlepy.org/docs/dev/) on [Tornado](http://www.tornadoweb.org/en/stable/). :book:

#### Tools

* [razi](https://github.com/rvianello/razi) - Cheminformatic extension for the SQLAlchemy database. :book:
* [ChemicaLite](https://github.com/rvianello/chemicalite) - SQLite cartridge depended on RDKit and SQLite. Stop maintained. :book:
* [PosgreSQL Cartridge](https://rdkit.readthedocs.io/en/latest/Cartridge.html) - RDKit ofiicial database cartridge. :book:

<a id="docking"></a>

### Docking

* [AutoDock Vina](http://vina.scripps.edu/) - Molecular docking and virtual screening. :book:
* [smina](https://sourceforge.net/projects/smina/) - Customized [AutoDock Vina](http://vina.scripps.edu/) to better support scoring function development and high-performance energy minimization. :book:
* [qvina](https://github.com/QVina/qvina) - Quick Vina 2 is a fast and accurate molecular docking tool, attained at accurately accelerating AutoDock Vina. QVINA: [DOI:10.1093/bioinformatics/btv082](https://doi.org/10.1093/bioinformatics/btv082), QVINA-W: [DOI:10.1038/s41598-017-15571-7](https://doi.org/10.1038/s41598-017-15571-7) :book:
* [rDock](http://rdock.sourceforge.net/) - rDock is a fast and versatile Open Source docking program that can be used to dock small molecules against proteins and nucleic acids. It is designed for High Throughput Virtual Screening (HTVS) campaigns and Binding Mode prediction studies. :book:
* [UCSF DOCK](http://dock.compbio.ucsf.edu/Overview_of_DOCK/index.htm) - DOCK 3.6 and DOC K6 are docking sofware designed by UCSF. Academic free. :credit_card:
* [Glide](https://www.schrodinger.com/Glide/) - A complete solution for ligand-receptor docking. Glide offers the full range of speed vs. accuracy options, from the HTVS, SP to XP. :moneybag:


<a id="md"></a>

### Molecular Dynamics

#### MD Engine

* [Gromacs](http://www.gromacs.org/) - Molecular dynamics package mainly designed for simulations of proteins, lipids and nucleic acids. :book:
* [OpenMM](http://openmm.org/) - High performance toolkit for molecular simulation including extensive language bindings for Python, C, C++, and even Fortran. :book:
* [NAMD](http://www.ks.uiuc.edu/Research/namd/) - NAMD is a parallel molecular dynamics code designed for high-performance simulation of large biomolecular systems.  :book:
* [Amber](http://ambermd.org/) - Amber is a suite of biomolecular simulation programs.  :moneybag:
* [Desmond](https://www.deshawresearch.com/resources_desmond.html) - Desmond is a software package developed at D. E. Shaw Research to perform high-speed molecular dynamics simulations of biological systems on conventional commodity clusters, general-purpose supercomputers, and GPUs. Non-commercial free. :moneybag:

#### Tools

* [CHARMM-GUI](http://charmm-gui.org/) - CHARMM-GUI has proven to be an ideal web-based platform to interactively build complex systems and prepare their inputs with well-established and reproducible simulation protocols for state-of-the-art biomolecular simulations using widely used simulation packages such as CHARMM, NAMD, GROMACS, AMBER, GENESIS, LAMMPS, Desmond, and OpenMM.  :book:
* [MDTraj](https://github.com/mdtraj/mdtraj) - Analysis of molecular dynamics trajectories. :book:
* [MSMBuilder](https://github.com/msmbuilder/msmbuilder) - MSMBuilder is a python package which implements a series of statistical models for high-dimensional time-series. It is particularly focused on the analysis of atomistic simulations of biomolecular dynamics. For example, MSMBuilder has been used to model protein folding and conformational change from molecular dynamics (MD) simulations.  :book:
* [PyEMMA](https://github.com/markovmodel/PyEMMA)PyEMMA (EMMA = Emma's Markov Model Algorithms) is an open source Python/C package for analysis of extensive molecular dynamics simulations and build Makov State Model. :book:

<a id="protein"></a>

### Protein
* [PDBFixer](https://github.com/pandegroup/pdbfixer) - PDBFixer is an easy to use application for fixing problems in Protein Data Bank files in preparation for simulating them. :book:
* [Biopython](https://biopython.org/) - Biopython is a set of freely available tools for biological computation written in Python by an international team of developers. :book:
* [ODDT](https://github.com/oddt/oddt) - Open Drug Discovery Toolkit (ODDT) is modular and comprehensive toolkit for use in cheminformatics, molecular modeling etc. ODDT is written in Python, and make extensive use of Numpy/Scipy. :book:
* [ProDy](http://prody.csb.pitt.edu) - ProDy is a free and open-source Python package for protein structural dynamics analysis, especially Elastic Network Model. :book:
* [Bio3D](thegrantlab.org/bio3d/index.php) - Bio3D is an R package containing utilities for the analysis of protein structure, sequence and trajectory data. :book:

## Journals

* [Journal of Cheminformatics](https://jcheminf.biomedcentral.com/)
* [Journal of Chemical Information and Modeling (ACS Publications)](https://pubs.acs.org/journal/jcisd8)

## Resources

### Courses

* [Learncheminformatics.com](http://learncheminformatics.com/) - "Cheminformatics: Navigating the world of chemical data" courese at Indiana University.
* [cheminfoeducation](https://www.youtube.com/user/cheminfoeducation/videos) - A YouTube channel for cheminformatics education.
* [CH485-Artificial-Intelligence-and-Chemistry](https://github.com/Birdlet/CH485---Artificial-Intelligence-and-Chemistry) - Course material for Artificial Intelligence and Chemistry of Korea Advanced Institute of Science and Technology (KAIST).
[CHEM430-SIMULATION IN CHEMISTRY & BIOCHEMISTRY](https://dasher.wustl.edu/chem430/) - IMULATION IN CHEMISTRY & BIOCHEMISTRY by Department of Chemistry at Washington University.
[CHEM478-MOLECULAR MODELING](https://dasher.wustl.edu/chem430/) - MOLECULAR MODELING by Department of Chemistry at Washington University.

### Blogs

* [Open Source Molecular Modeling](https://opensourcemolecularmodeling.github.io/README.html) - Updateable catalog of open source molecular modeling software.
* [PubChem Blog](https://pubchemblog.ncbi.nlm.nih.gov/) - News, updates and tutorials about [PubChem](https://pubchem.ncbi.nlm.nih.gov/).
* [The ChEMBL-og blog](http://chembl.blogspot.tw/) - Stories and news from Computational Chemical Biology Group at [EMBL-EBI](https://www.ebi.ac.uk/).
* [ChEMBL blog](http://chembl.github.io/) - ChEMBL on GitHub.
* [SteinBlog](http://www.steinbeck-molecular.de/steinblog/) - Blog of [Christoph Steinbeck](http://www.steinbeck-molecular.de/steinblog/index.php/about/), who is the head of cheminformatics and metabolism at the EMBL-EBI.
* [Practical Cheminformatics](http://practicalcheminformatics.blogspot.com/) - Blog with in-depth examples of practical application of cheminformatics.
* [So much to do, so little time - Trying to squeeze sense out of chemical data](http://blog.rguha.net/) - Bolg of [Rajarshi Guha](http://blog.rguha.net/?page_id=8), who is a research scientist at NIH Center for Advancing Translational Science.
  * Some old blogs [1](https://rguha.wordpress.com/) [2](http://www.rguha.net/index.html).
* [Noel O'Blog](http://baoilleach.blogspot.tw/) - Blog of [Noel O'Boyle](https://www.redbrick.dcu.ie/~noel/), who is a Senior Software Engineer at NextMove Software.
* [chem-bla-ics](http://chem-bla-ics.blogspot.tw/) - Blog of [Egon Willighagen](http://egonw.github.io/), who is an assistant professor at Maastricht University.
* [Asad's Blog](https://chembioinfo.com/) - Bolg of Syed Asad Rahman, who is a research scientist in the [Thornton group](http://www.ebi.ac.uk/research/thornton) at EMBL-EBI.
* [steeveslab-blog](http://asteeves.github.io/) - Some examples using [RDKit](http://www.rdkit.org/).
* [Macs in Chemistry](http://www.macinchem.org/) - Provide a resource for chemists using Apple Macintosh computers.
* [Is Life Worth Living](https://iwatobipen.wordpress.com) - iwatobipen is a medicinal chemist in mid-size pharmaceutical company in Japan. He loves cheminformatics and posted lots of chemoinformatical blogs

### Books

* [Computational Approaches in Cheminformatics and Bioinformatics](https://books.google.com/books/about/Computational_Approaches_in_Cheminformat.html?id=bLqV4rYQoYsC) -  Include insights from public (NIH), academic, and industrial sources at the same time.
* [Chemoinformatics for Drug Discovery](https://onlinelibrary.wiley.com/doi/book/10.1002/9781118742785) - Materials about how to use Chemoinformatics strategies to improve drug discovery results.
* [Handbook of Chemoinformatics: From Data to Knowledge in 4 Volumes](https://onlinelibrary.wiley.com/doi/book/10.1002/9783527618279) - Many detials for cheminformatics descriptors and algorithm. RDkit cited this book alot.

<a id="see-also"></a>
## See Also

* [deeplearning-biology](https://github.com/hussius/deeplearning-biology#chemoinformatics-and-drug-discovery-) - Chemoinformatics and drug discovery section in deeplearning-biology repo.
* [awesome-python-chemistry](https://github.com/lmmentel/awesome-python-chemistry) - Another list focuses on Python stuff related to Chemistry.

## License

[![CC0](http://mirrors.creativecommons.org/presskit/buttons/88x31/svg/cc-zero.svg)](https://creativecommons.org/publicdomain/zero/1.0/)
