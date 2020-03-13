% Runner Class that controls the larger Analysis Class.
%
% Author: Joshua Rees, joshua.rees@bristol.ac.uk
% Affiliation: BrisSynBio, Life Sciences, University of Bristol
% Last Updated: 07 / 03 / 2017
% Modified from Original Work: https://github.com/CovertLab/WholeCell/blob/master/src/%2Bedu/%2Bstanford/%2Bcovert/%2Bcell/%2Bsim/%2Banalysis/SingleGeneDeletions.m
% Original Author: Jonathan Karr

function runWholeCellAnalysis()
    %setWarnings();
    %setPath();

    % change to MY CLASS
    %import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;

    %Declarations
    % NO: baseDir
    % 'C:', 'Users', 'jr0904', 'Google Drive', 'COMPUTER', 'WholeCell-master', 'WholeCell-master'
    % .../projects/TestRun/
    % .../projects/TestRun/AnalysisOutput
    % .../projects/WildType/1
    % .../projects/TestRun/RawData/1 - 5
    % CHANGE everytime as hardcoded, and CREATE folders/data

    %fileName - set in load, reminder to change, and what is>here
    projectDir = fullfile('C:', 'Users', 'jr0904', 'Google Drive', 'COMPUTER', 'WholeCell-master', 'WholeCell-master',...
        'projects', 'TestRun');
    outputDir = [projectDir filesep 'AnalysisOutput'];
    wildtypeDir = [projectDir filesep 'WildType' filesep '1'];
    rawdataDir = [projectDir filesep 'RawData'];
    simulationNumber = 1;
    numberofSimulations = 5;
    %looped below - simulationDir = [rawdataDir filesep num2str(simulationNumber)];
    wildtypeFileName = 'projects/TestRun/wildtype-fitted%s.mat';
    %looped below - simulationFileName = sprintf('projects/TestRun/simulation-%d-fitted%s.mat', simulationNumber);

    % ADD CHECKS to see if files already exist ...
    % ADD CHECKS to see if metadata.mat exists ... > how check mutants
    % otherwise?

    %Loading of WC output data, processing and storage
    directoryOrWID = wildtypeDir;
    processRawData(directoryOrWID, wildtypeFileName);

    %Loading of Simulations output data, processing and storage
    for sims = simulationNumber:numberofSimulations
        simulationDir = [rawdataDir filesep num2str(simulationNumber)];
        simulationFileName = sprintf('projects/TestRun/simulation-%d-fitted%s.mat', simulationNumber);

        directoryOrWID = simulationDir;
        processRawData(directoryOrWID, simulationFileName);

        simulationNumber = simulationNumber + 1;
    end

    %Recreating figure 6 from Karr Paper
    % change to MY OWN CLASS
    % check args, use cacheddata, geneIDX etc

    %SingleGeneDeletions. 'Classify Inidividual Deletion Strains'
    %SingleGeneDeletions.plotOverview(...
    %[outDir filesep 'projectname.pdf'], ...
    %[outDir filesep 'projectname.xls']);
end

%DOCUMENTATION
    %GETTING STARTED WITH WholeCell Modelling
        % DOWNLOAD WholeCell from https://simtk.org/projects/wholecell
        % INSTALL according to the tutorial http://web.archive.org/web/20160218134737/http://www.wholecell.org/wiki/index.php?title=Tutorial
        % BASE DIR is suggested as /WholeCell-master/output/runSimulation
        % EXPECTED folder structure /runSimulation/yyyy_mm_dd_HH_MM_SS/1 ...(to N) 
        % RUN simulations locally in Matlab using tutorial Box1 and Box3
        % PLOT Box 6 (changing file names or assigning cached data>'log')
        % PLOT Box 7 using the %import classes, %load data, %plot sections...
        % ...plus declare simBatchDir (see below) and simIdx = 1; in your script
        % RUN simulations on a cluster using our slurm scripts
        % ANALYSE local data or online data downloaded ...
        % ...to EXPECTED folder structure, using this script.
