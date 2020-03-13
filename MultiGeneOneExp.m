% MultigeneOneExp class, controlled by Runner File
%
% Author: Joshua Rees, joshua.rees@bristol.ac.uk
% Affiliation: BrisSynBio, Life Sciences, University of Bristol
% Last Updated: 07 / 03 / 2017
% Modified from Original Work: https://github.com/CovertLab/WholeCell/blob/master/src/%2Bedu/%2Bstanford/%2Bcovert/%2Bcell/%2Bsim/%2Banalysis/SingleGeneDeletions.m
% Original Author: Jonathan Karr

%CHANGE SingleGeneDeletions all at once
classdef SingleGeneDeletions
    %Properties taken from Karr's SingleGeneDeletions.m
    properties (Constant = true)
        %A. Non-essential (NON_ESSENTIAL)
        %B. Semi-essential
        %   1. Slow growing (SLOW_GROWING)
        %C. Essential
        %   1. No growth
        %      a. Macromolecule maintenance (NON_GROWING)
        %      b. No macromolecule maintenance, toxic metabolite levels (DECOMPOSING)
        %   2. Decaying growth
        %      a. Non-RNA, Non-protein synthesizing (DECAYING_GROWTH_NON_RNA_PROTEIN)
        %         i.  Non-transcribing
        %         ii. Non-maturing
        %      b. Non-protein synthesizing (DECAYING_GROWTH_NON_PROTEIN)
        %   3. Non-dividing
        %      a. Non-replicative (NON_REPLICATIVE)
        %      b. Non-fissive (NON_FISSIVE)
        %   4. Unable to sustain division over many generations (NON_PERPETUATING)
        %   5. No terminal organelle (NO_TERMINAL_ORGANELLE)
        %   6. Toxin accumulation (TOXIN_ACCUMULATION)
        DELETION_STRAIN_CLASSES = {
            'Non-essential'             'Doesn''t meet criterion of any of the other categories'
            'Decomposing'               'No growth, no energy production, no macromolecule maintenance, toxic metabolite levels'
            'Non-growing'               'No growth, energy production, macromolecule maintenance'
            'Non-RNA synthesizing'      'Derivative of growth is negative, time to mass doubling is long, time to end of replication is long, cell cycle is elongated'
            'Non-protein synthesizing'  'Derivative of growth is negative, time to mass doubling is long, time to end of replication is long, cell cycle is elongated'
            'Non-replicative'           'End ploidy = 1, high dNTPs, no simulation finishes replication (rep end time = NaN)'
            'Non-fissive'               'End ploidy = 2, but pinchedDiamter = diameter'
            'Slow growing'              'Grows and divides, but slower than wild-type'
            'Toxin accumulation'        'Accumulates toxic concentrations of metabolites'
            'Non-perpetuating'          'Gene is essential, but cells will survive for several generations until gene products are sufficiently diluted among progeny'
            'No terminal organelle'     'Lack terminal organelle'
            'Unobserved'                'Simulation not yet run'
            'Wild-type'                 'Wild-type simulation'
            };
        NON_ESSENTIAL                   = 1
        DECOMPOSING                     = 2
        NON_GROWING                     = 3
        DECAYING_GROWTH_NON_RNA_PROTEIN = 4
        DECAYING_GROWTH_NON_PROTEIN     = 5
        NON_REPLICATIVE                 = 6
        NON_FISSIVE                     = 7
        SLOW_GROWING                    = 8
        TOXIN_ACCUMULATION              = 9
        NON_PERPETUATING                = 10
        NO_TERMINAL_ORGANELLE           = 11
        UNOBSERVED                      = 12
        WILD_TYPE                       = 13
        
        COLORS = [
            0 0 1 %blue
            1 0 0 %red
            0 1 0 %green
            0 1 1 %cyan
            1 0.5 0 %orange
            0.5 0 0.5 %purple
            1 0.84 0 %gold
            0.25 0 0.5 %indigo
            0 1 1 %aqua
            0.5 1 0.8314 %aquamarine
            0.329 0.525 0.043 %dark gold
            1 0.753 0.795 %pink
            0.25 0.25 0.25 %dark grey
            ];
        
        nonPerpetuatingEssentialGenes = {
            'MG_0001'
            'MG_012'
            'MG_019'
            'MG_048'
            'MG_072'
            'MG_106'
            'MG_109'
            'MG_110'
            'MG_139'
            'MG_143'
            'MG_170'
            'MG_172'
            'MG_184'
            'MG_210'
            'MG_277'
            'MG_297'
            'MG_305'
            'MG_329'
            'MG_335'
            'MG_384'
            'MG_387'
            'MG_392'
            'MG_393'
            'MG_425'
            'MG_442'
            'MG_464'
            'MG_476'
            };
    end
    methods (Static)
        
        %Classify Deletion Strains
        % Recieves a declared mutant mat file, compares it to WT data and 
        % determines the class of the mutation.
        % change basedir, deletiondata, wtdata
        function [class, classLabel] = classifyDeletionStrain(geneID, wtData)
        % geneID, wtData args passed in, class, classlabel returned out
            
        %import classes
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %declare variables / load data from args
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            deletionData = load([baseDir geneID '.mat']);
            if nargin < 2
                wtData = load([baseDir 'WT.mat']);
            end
            nWtSims = numel(wtData.simEndTimes);
            nDeletionSims = numel(deletionData.simEndTimes);
            cellCycleLength = min([wtData.simEndTimes; deletionData.simEndTimes]);
            
            %determining the class of mutation, which is returned
            wtStartTimes = 2 * ones(size(wtData.simEndTimes));
            wtTauTimes = zeros(size(wtData.growth, 2), 1);
            wtFinTimes = wtData.simEndTimes + 1;
            deletionStartTimes = 101 * ones(size(deletionData.simEndTimes));
            deletionTauTimes = zeros(size(deletionData.growth, 2), 1);
            deletionFinTimes = deletionData.simEndTimes + 1;
            for i = 1:nWtSims
                wtTauTimes(i) = find(~isnan(wtData.growth(1:cellCycleLength, i)), 1, 'last');
            end
            for i = 1:nDeletionSims
                deletionTauTimes(i) = find(~isnan(deletionData.growth(1:cellCycleLength, i)), 1, 'last');
            end
            
            wtStartInds = sub2ind(size(wtData.growth), wtStartTimes', 1:nWtSims);
            wtFinInds = sub2ind(size(wtData.growth), wtFinTimes', 1:nWtSims); %#ok<NASGU>
            wtTauInds = sub2ind(size(wtData.growth), wtTauTimes', 1:nWtSims);
            deletionStartInds = sub2ind(size(deletionData.growth), deletionStartTimes', 1:nDeletionSims);
            deletionFinInds = sub2ind(size(deletionData.growth), deletionFinTimes', 1:nDeletionSims);
            deletionTauInds = sub2ind(size(deletionData.growth), deletionTauTimes', 1:nDeletionSims);
            
            wtTauGrowths = wtData.growth(wtTauInds);
            deletionFinGrowths = deletionData.growth(deletionFinInds);
            deletionTauGrowths = deletionData.growth(deletionTauInds);
            deletionInitGrowths = deletionData.growth(deletionStartInds);
            
            if ...
                    median(deletionInitGrowths) < 0.1 * median(wtData.growth(wtStartInds)) && (...
                    median(deletionData.damagedProteins(deletionFinInds)) > max(100, 100 * median(wtData.damagedProteins(wtTauInds))) || ...
                    median(deletionData.damagedRnas(deletionFinInds))     > max(100, 100 * median(wtData.damagedRnas(wtTauInds)))     || ...
                    median(deletionData.aminoAcids(deletionFinInds))      > max(100, 100 * median(wtData.aminoAcids(wtTauInds)))      || ...
                    median(deletionData.ntps(deletionFinInds))            > max(100, 100 * median(wtData.ntps(wtTauInds)))               ...
                    )
                class = SingleGeneDeletions.DECOMPOSING;
            elseif median(deletionInitGrowths) < 0.1 * median(wtData.growth(wtStartInds))
                class = SingleGeneDeletions.NON_GROWING;
            elseif all(deletionFinGrowths < deletionInitGrowths) && all(deletionData.rnas(deletionTauInds) < deletionData.rnas(deletionStartInds))
                class = SingleGeneDeletions.DECAYING_GROWTH_NON_RNA_PROTEIN;
            elseif all(deletionFinGrowths < deletionInitGrowths) && all(deletionData.prots(deletionTauInds) < deletionData.prots(deletionStartInds))
                class = SingleGeneDeletions.DECAYING_GROWTH_NON_PROTEIN;
            elseif ...
                    median(deletionData.damagedProteins(deletionFinInds)) > max(100, 100 * median(wtData.damagedProteins(wtTauInds))) || ...
                    median(deletionData.damagedRnas(deletionFinInds))     > max(100, 100 * median(wtData.damagedRnas(wtTauInds)))     || ...
                    median(deletionData.aminoAcids(deletionFinInds))      > max(100, 100 * median(wtData.aminoAcids(wtTauInds)))      || ...
                    median(deletionData.ntps(deletionFinInds))            > max(100, 100 * median(wtData.ntps(wtTauInds)))            || ...
                    median(deletionData.antibiotics(deletionFinInds))     > max(100, 100 * median(wtData.antibiotics(wtTauInds)))
                class = SingleGeneDeletions.TOXIN_ACCUMULATION;
            elseif ~any(deletionData.ploidy(:) > 1)
                class = SingleGeneDeletions.NON_REPLICATIVE;
            elseif median(deletionTauGrowths) < 0.80 * median(wtTauGrowths) && median(deletionTauGrowths) > median(deletionInitGrowths)
                class = SingleGeneDeletions.SLOW_GROWING;
            elseif ~any(deletionData.cytokinesisDuration(:))
                class = SingleGeneDeletions.NON_FISSIVE;
            elseif ~any(deletionData.isBacteriumAdherent(2, :))
                class = SingleGeneDeletions.NO_TERMINAL_ORGANELLE;
            elseif ismember(geneID, SingleGeneDeletions.nonPerpetuatingEssentialGenes)
                class = SingleGeneDeletions.NON_PERPETUATING;
            else
                class = SingleGeneDeletions.NON_ESSENTIAL;
            end
            
            % determining the class label, which is returned
            classLabel = SingleGeneDeletions.DELETION_STRAIN_CLASSES{class};
        end
        
        
        %Classify Individual Deletion Strains
        % Creates geneclasses.mat by calling classifydeletionstrains and
        % cycling mutant mat files.
        % change name to classifyMultipleDeletionStrains
        % use of geneIDXs
        % change basedir, sim, pass in project filename example (>g.WholeCellModelIDS),  
        % change WT.mat, MG*.mat
        function geneClasses = classifyIndividualDeletionStrains(useCachedData)
        % useCachedData (absence or presence) input, geneClasses.mat stored
            %import classes
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %declare variables / load data from args
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            sim = CachedSimulationObjectUtil.load();
            g = sim.gene;
            wtData = load([baseDir 'WT.mat']);
            tmp = dir([baseDir filesep 'MG*.mat']);
            
            %set up loop to push mutant files to get classfied and saved in gene classes 
            [~, geneticKnockouts] = ismember(cellfun(@(str) str(1:end-4), {tmp.name}', 'UniformOutput', false), g.wholeCellModelIDs);
            geneIdxs = setdiff(unique(geneticKnockouts), 0);
            
            %check see if geneclasses exist and if it has any missing data
            if nargin >= 1 && useCachedData && exist([baseDir 'geneClasses.mat'], 'file')
                load([baseDir 'geneClasses.mat']);
                geneIdxs = intersect(geneIdxs, find(geneClasses == SingleGeneDeletions.UNOBSERVED)); %#ok<NODEF>
            else
                %creating empty geneclasses / preallocation?
                geneClasses = repmat(SingleGeneDeletions.UNOBSERVED, size(g.wholeCellModelIDs));
            end
            
            %replace this loop (remove geneIdxs)
            for i = 1:numel(geneIdxs)
                geneIdx = geneIdxs(i);
                geneID = g.wholeCellModelIDs{geneIdx};
                geneClasses(geneIdx) = SingleGeneDeletions.classifyDeletionStrain(geneID, wtData);
                clear deletionData;
            end
            
            try %#ok<TRYNC>
                save([baseDir 'geneClasses.mat'], 'geneClasses');
            end
        end
        
        %getgeneimplementations
        % tied into geneIDXs in CalcOverview. Couldn't tell you what it
        % does. May not need it.
        function [tfs, molecules, processes] = getGeneImplementations(sim)
            g = sim.gene;
            
            tfs = false(size(g.wholeCellModelIDs));
            molecules = repmat({{}}, size(g.wholeCellModelIDs));
            processes = repmat({{}}, size(g.wholeCellModelIDs));
            
            for i = 1:length(sim.processes)
                m = sim.processes{i};
                
                %stimuli
                stimuliGeneComposition = m.stimuliGeneComposition();
                for j = 1:length(sim.gene.wholeCellModelIDs)
                    for k = 1:length(m.stimuliWholeCellModelIDs)
                        if ~stimuliGeneComposition(j, k)
                            continue;
                        end
                        
                        tfs(j) = true;
                        molecules{j} = [molecules{j}; m.stimuliWholeCellModelIDs(k)];
                        processes{j} = [processes{j}; m.wholeCellModelID];
                    end
                end
                
                %substrates
                substrateGeneComposition = m.substrateGeneComposition();
                for j = 1:length(sim.gene.wholeCellModelIDs)
                    for k = 1:length(m.substrateWholeCellModelIDs)
                        if ~substrateGeneComposition(j, k)
                            continue;
                        end
                        
                        tfs(j) = true;
                        molecules{j} = [molecules{j}; m.substrateWholeCellModelIDs(k)];
                        processes{j} = [processes{j}; m.wholeCellModelID];
                    end
                end
                
                %enzymes
                enzymeGeneComposition = m.enzymeGeneComposition();
                for j = 1:length(sim.gene.wholeCellModelIDs)
                    for k = 1:length(m.enzymeWholeCellModelIDs)
                        if ~enzymeGeneComposition(j, k)
                            continue;
                        end
                        
                        tfs(j) = true;
                        molecules{j} = [molecules{j}; m.enzymeWholeCellModelIDs(k)];
                        processes{j} = [processes{j}; m.wholeCellModelID];
                    end
                end
            end
        end
        
        %calcOverviewData
        % creates gridData, gridMetaData (=Overview.mat) and Overview-....mat
        % called by plotOverviewData
        % In Future: debate over whether need - groups sims together around class of
        % mutation - may not be quite what want
        
        % biggest implementation of geneidx
        % change basedir to datadir, ...
        function [gridData, gridMetaData] = calcOverviewData(iJob, nJobs)
        %currently input args are unused, output gridData, gridMetaData
            
            %import classes
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %declarations
            if nargin < 2
                iJob = 1;
                nJobs = 1;
            end
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            load([baseDir 'geneClasses.mat']);
            sim = CachedSimulationObjectUtil.load();
            g = sim.gene;
            isGeneImplemented = SingleGeneDeletions.getGeneImplementations(sim);
            nCats = size(SingleGeneDeletions.DELETION_STRAIN_CLASSES, 1) - 1;
            nProps = 8;
            nTime = 50001;
            catLabels = {
                'WT'             SingleGeneDeletions.WILD_TYPE                        'growth'
                'Energy'         SingleGeneDeletions.DECOMPOSING                      'ntps'
                'Metabolic'      SingleGeneDeletions.NON_GROWING                      'growth'
                'RNA'            SingleGeneDeletions.DECAYING_GROWTH_NON_RNA_PROTEIN  'rnaWt'
                'Protein'        SingleGeneDeletions.DECAYING_GROWTH_NON_PROTEIN      'proteinWt'
                'Other'          SingleGeneDeletions.NON_PERPETUATING                 'proteinWt'
                'DNA'            SingleGeneDeletions.NON_REPLICATIVE                  'dnaWt'
                'Cytokinesis'    SingleGeneDeletions.NON_FISSIVE                      'growth'
                'Term Org'       SingleGeneDeletions.NO_TERMINAL_ORGANELLE            'terminalOrganelleWt'
                'Damaged'        SingleGeneDeletions.TOXIN_ACCUMULATION               'damagedProteins'
                'Slow'           SingleGeneDeletions.SLOW_GROWING                     'rnaWt'
                'Non-Ess'        SingleGeneDeletions.NON_ESSENTIAL                    'growth'
                };
            propLabels = {
                {'NTP'}                 'ntps'
                {'Growth (fg h^{-1})'}  'growth'
                {'Protein (fg)'}        'proteinWt'
                {'RNA (fg)'}            'rnaWt'
                {'DNA (fg)'}            'dnaWt'
                {'Septum (nm)'}         'pinchedDiameter'
                {'Term Org (fg)'}       'terminalOrganelleWt'
                {'Damaged' 'Prot'}      'damagedProteins'
                };
            
            %initialising gridData and gridMetaData
            gridData = NaN(nProps, nCats, nTime);
            gridMetaData = repmat(struct('nSimulations', [], 'nGenes', [], 'gene', [], 'simGroup', [], 'simIdx', []), nCats, 1);
            
            %filling the gridData via tmp / tmp2
            for i = iJob:nJobs:nCats
                if exist([baseDir 'overview-' catLabels{i, 1} '.mat'], 'file')
                    load([baseDir 'overview-' catLabels{i, 1} '.mat']);
                else
                    %get data for all simulations in category
                    % propLabels WT
                    if catLabels{i, 2} == SingleGeneDeletions.WILD_TYPE
                        tmp2 = load([baseDir 'WT.mat'], 'metaData', 'simEndTimes', 'ntps', 'growth', 'rnaWt', 'proteinWt', 'dnaWt', 'pinchedDiameter', 'terminalOrganelleWt', 'damagedProteins');
                        tmp = tmp2.growth(sub2ind(size(tmp2.growth), tmp2.simEndTimes, (1:numel(tmp2.simEndTimes))'));
                        for k = 1:nProps
                            for l = 1:numel(tmp2.simEndTimes)
                                if isnan(tmp(l)) || ~isfield(tmp2, propLabels{k, 2}) || any(isnan(tmp2.(propLabels{k, 2})(2:tmp2.simEndTimes(l), l)))
                                    tmp(l) = NaN;
                                end
                            end
                        end
                        [~, idx] = min(abs(tmp - nanmedian(tmp)));
                        
                        tmpMetaData = struct(...
                            'nSimulations', numel(tmp), ...
                            'nGenes', [], ...
                            'gene', [], ...
                            'simGroup', tmp2.metaData(idx).simGroup, ...
                            'simIdx', tmp2.metaData(idx).simIdx);
                    % catLabels mutant
                    else
                        geneIdxs = find(geneClasses == catLabels{i, 2} & isGeneImplemented);
                        
                        tmp = zeros(0, 1);
                        genes = zeros(0, 1);
                        simIdxs = zeros(0, 1);
                        for j = 1:numel(geneIdxs)
                            tmp2 = load([baseDir g.wholeCellModelIDs{geneIdxs(j)} '.mat'], 'simEndTimes', catLabels{i, 3});
                            if ~isfield(tmp2, catLabels{i, 3})
                                continue;
                            end
                            tmp = [tmp; tmp2.(catLabels{i, 3})(sub2ind(size(tmp2.(catLabels{i, 3})), tmp2.simEndTimes, (1:numel(tmp2.simEndTimes))'))];
                            genes = [genes; repmat(geneIdxs(j), numel(tmp2.simEndTimes), 1)];
                            simIdxs = [simIdxs; (1:numel(tmp2.simEndTimes))'];
                        end
                        
                        % propLabels mutant
                        while true
                            switch catLabels{i, 2}
                                case SingleGeneDeletions.NO_TERMINAL_ORGANELLE
                                    [~, idx] = min(tmp);                                
                                otherwise
                                    [~, idx] = min(abs(tmp - nanmedian(tmp)));
                            end
                            tmp2 = load([baseDir g.wholeCellModelIDs{genes(idx)} '.mat'], 'metaData', 'simEndTimes', propLabels{:, 2});
                            
                            if ~all(ismember(propLabels(:, 2), fieldnames(tmp2)))
                                tmp(idx) = NaN;
                            end
                            for k = 1:nProps
                                if ~isfield(tmp2, propLabels{k, 2}) || any(isnan(tmp2.(propLabels{k, 2})(2:tmp2.simEndTimes(simIdxs(idx)), simIdxs(idx))))
                                    tmp(idx) = NaN;
                                end
                            end
                            
                            if catLabels{i, 2} == SingleGeneDeletions.NON_FISSIVE && (...
                                    isnan(range(tmp2.dnaWt(:, simIdxs(idx)))) || ~range(tmp2.dnaWt(:, simIdxs(idx))) || ...
                                    isnan(tmp2.rnaWt(end-4*3600, simIdxs(idx))) || tmp2.rnaWt(end, simIdxs(idx)) < tmp2.rnaWt(end-4*3600, simIdxs(idx)))
                                tmp(idx) = NaN;
                            end
                            
                            if ~isnan(tmp(idx))
                                break;
                            end
                        end
                        
                        tmpMetaData = struct(...
                            'nSimulations', numel(tmp), ...
                            'nGenes', numel(geneIdxs), ...
                            'gene', g.wholeCellModelIDs{genes(idx)}, ...
                            'simGroup', tmp2.metaData(simIdxs(idx)).simGroup, ...
                            'simIdx', tmp2.metaData(simIdxs(idx)).simIdx);
                        
                        idx = simIdxs(idx);
                    end
                    
                    %get data for representative cell
                    tmpData = struct;
                    for k = 1:nProps
                        tmpData.(propLabels{k, 2}) = tmp2.(propLabels{k, 2})(:, idx);
                    end
                    
                    %save
                    save([baseDir 'overview-' catLabels{i, 1} '.mat'], 'tmpMetaData', 'tmpData');
                    
                    %cleanup
                    clear tmp2;
                end
                
                %gridMetaData
                gridMetaData(i) = tmpMetaData;
                %gridData
                for j = 1:nProps
                    gridData(j, i, 1:numel(tmpData.(propLabels{j, 2}))) = tmpData.(propLabels{j, 2});
                end
            end
        end
        
        
        %plotOverviewData
        % main plotting function, producing figure 6 graph
        % change inputs and outputs
        % remove figures don't want
        % heavily customised. breakdown and customise after first run.
        function plotOverview(outFileName_Figure, outFileName_Table, ...
                sim, geneClasses, gridData, gridMetaData)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            if nargin < 1 || isempty(outFileName_Figure)
                outFileName_Figure = [baseDir filesep 'overview.pdf'];
            end
            if nargin < 2 || isempty(outFileName_Table)
                outFileName_Table = [baseDir 'overview.xls'];
            end
            
            %% initialize
            figW = 17.4;
            figH = 7.3;
            
            [~, figHandle] = PlotUtil.newAxesHandle();
            clf(figHandle);
            set(figHandle, 'PaperUnits', 'centimeters');
            set(figHandle, 'PaperSize', [figW figH]); %cell two column
            set(figHandle, 'PaperPositionMode', 'manual');
            set(figHandle, 'PaperPosition', [0 0 get(figHandle, 'PaperSize')]);
            
            %% get data
            if nargin < 3
                sim = CachedSimulationObjectUtil.load();
            end
            g = sim.gene;
            isGeneImplemented = SingleGeneDeletions.getGeneImplementations(sim);
            expEss = ismember(g.essential, {'Y', ''});
            if nargin < 4
                load([baseDir 'geneClasses.mat']);
            end
            modelEss = geneClasses ~= SingleGeneDeletions.NON_ESSENTIAL;
            
            showGeneSymbols = true;
            geneSymbols = struct;
            geneSymbols.MG006 = 'tmk';
            geneSymbols.MG022 = 'rpoE';
            geneSymbols.MG113 = 'asnS';
            geneSymbols.MG048 = 'ffh';
            geneSymbols.MG001 = 'dnaN';
            geneSymbols.MG204 = 'parC';
            geneSymbols.MG084 = 'tilS';
            
            %% sizes
            fontSizeSubfigLabel = 8;
            fontSizeLarge = 7;
            fontSizeMed = 6;
            fontSizeSmall = 5;
            
            leftColX = 0.05;
            leftColLabelX = 0.051;
            leftColW = 0.15;
            leftColH1 = .23;
            leftColH2 = 0.23;
            leftColH3 = leftColH2;
            leftColY2 = 0.375;
            leftColY3 = 0.0622;
            
            %% part A -- model vs. observed gene essentiality
            %PlotUtil.labelSubFigure('A', [leftColLabelX 1.01 -0.0286 -0.0690], figHandle);
            %SingleGeneDeletions.plotGeneDeletionPrediction(sim, geneClasses, [leftColX 0.71 leftColW leftColH1]);
            
            %% part B -- deletion strain classes spark line matrix
            catLabels = {
                'WT'             SingleGeneDeletions.WILD_TYPE
                'Energy'         SingleGeneDeletions.DECOMPOSING
                'Metabolic'      SingleGeneDeletions.NON_GROWING
                'RNA'            SingleGeneDeletions.DECAYING_GROWTH_NON_RNA_PROTEIN
                'Protein'        SingleGeneDeletions.DECAYING_GROWTH_NON_PROTEIN
                'Other'          SingleGeneDeletions.NON_PERPETUATING
                'DNA'            SingleGeneDeletions.NON_REPLICATIVE
                'Cytokinesis'    SingleGeneDeletions.NON_FISSIVE
                'Term Org'       SingleGeneDeletions.NO_TERMINAL_ORGANELLE
                'Damaged'        SingleGeneDeletions.TOXIN_ACCUMULATION
                'Quasi-Ess'      SingleGeneDeletions.SLOW_GROWING
                'Non-Ess'        SingleGeneDeletions.NON_ESSENTIAL
                };
            propLabels = {
                {'NTP'}                 'ntps'
                {'Growth (fg h^{-1})'}  'growth'
                {'Protein (fg)'}        'proteinWt'
                {'RNA (fg)'}            'rnaWt'
                {'DNA (fg)'}            'dnaWt'
                {'Septum (nm)'}         'pinchedDiameter'
                {'Term Org (fg)'}       'terminalOrganelleWt'
                {'Damaged' 'Prot'}      'damagedProteins'
                };
            highlightPlot = false(size(propLabels, 1), size(catLabels, 1));
            %update these, tried to preempt errors
            %highlightPlot(1:8, 2) = true;%energy
            highlightPlot(2:7, 3) = true;%metabolic
            highlightPlot([1:4 6:7], 4) = true;%RNA
            highlightPlot([2 3 6:7], 5) = true;%Protein
            highlightPlot(7, 6) = true;%other synthetic
            highlightPlot(5:6, 7) = true;%DNA
            highlightPlot(6, 8) = true;%cytokinesis
            %highlightPlot(8, 10) = true;%damaged
            %highlightPlot([2 3 7], 11) = true;%Slow
            
            groupCats = (1:size(catLabels, 1))';
            groupCats([2 9 10 12]) = [3 6 11 0];
            showProps = (1:size(propLabels, 1))';
            showProps([1 7 8]) = false;
            
            catIdxs = find(groupCats == (1:size(catLabels, 1))');
            propIdxs = find(showProps);
            
            nCats = numel(catIdxs);
            nProps = numel(propIdxs);
            
            time = (0:50000)' / 3600;
            
            if nargin >= 8
            elseif exist([baseDir 'overview.mat'], 'file')
                load([baseDir 'overview.mat']);
            else
                [gridData, gridMetaData] = SingleGeneDeletions.calcOverviewData();
                save([baseDir 'overview.mat'], 'gridData', 'gridMetaData');
            end
            %WHY save the damn thing if you don't load it again and just
            %give gridData a value???
            gridData(:, :, 1) = NaN;
            
            %layout
            x = leftColX + leftColW + 0.13;
            y = 0.845;
            W = 1 - x - 0.005;
            H = 0.773;
            w = W / nCats;
            h = H / nProps;
            y1 = y + 0.12;
            y2 = y + 0.08;
            y3 = y + 0.06;
            yMargin = 0.13 / 2;
            xMargin = (h * figH) / (w * figW) * yMargin;
            xSpan = 1 - 2 * xMargin;
            ySpan = 1 - 2 * yMargin;
            highlightColor = [255 230 230] / 255; %light red
            nonhighlightColor = [1 1 1]; %white
            
            %label
            PlotUtil.labelSubFigure('B', [x-0.0650 1.01 -0.0286 -0.0690], figHandle, fontSizeSubfigLabel);
            
            %data
            xlims = [0 time(end)] + [-1 1] * (1/xSpan - 1)/2 * time(end);
            for i = 1:nProps
                for j = 1:nCats
                    axesHandle = subplot('Position', [x+(j-1)*w  y-(i)*h  w  h]);
                    
                    switch propLabels{propIdxs(i), 2}
                        case 'ntps'
                            yticks = [0 1e6];
                            ytickLabels = {num2str(yticks(1))  '10^6'};
                        case 'growth'
                            yticks = [0 2.5];
                            ytickLabels = {num2str(yticks(1))  num2str(yticks(2))};
                        case 'rnaWt'
                            yticks = [0 0.4] * 1e-15;
                            ytickLabels = {num2str(yticks(1)*1e15)  num2str(yticks(2)*1e15)};
                        case 'proteinWt'
                            yticks = [2 7] * 1e-15;
                            ytickLabels = {num2str(yticks(1)*1e15)  num2str(yticks(2)*1e15)};
                        case 'dnaWt'
                            yticks = [0.6 1.2] * 1e-15;
                            ytickLabels = {num2str(yticks(1)*1e15)  num2str(yticks(2)*1e15)};
                        case 'pinchedDiameter'
                            yticks = [0 250e-9];
                            ytickLabels = {num2str(yticks(1))  num2str(yticks(2) * 1e9)};
                        case 'terminalOrganelleWt'
                            yticks = [0 0.1] * 1e-15;
                            ytickLabels = {num2str(yticks(1) * 1e15)  num2str(yticks(2) * 1e15)};
                        case 'damagedProteins'
                            yticks = [0 6000];
                            ytickLabels = {num2str(yticks(1))  num2str(yticks(2))};
                        otherwise
                            throw(MException('SingleGeneDeletions:error', 'undefined property %s', propLabels{propIdxs(i), 2}))
                    end
                    ylims0 = [min(yticks(1), min(min(gridData(propIdxs(i), :, :))))  max(yticks(2), max(max(gridData(propIdxs(i), :, :))))];
                    ylims = ylims0 + [-1 1] * (1/ySpan - 1)/2 * (ylims0(2) - ylims0(1));
                    
                    if highlightPlot(propIdxs(i), catIdxs(j))
                        faceColor = highlightColor;
                    else
                        faceColor = nonhighlightColor;
                    end
                    if ~all(faceColor == 1)
                        hold(axesHandle, 'on');
                        rectangle('Parent', axesHandle, 'Position', [xlims(1) ylims(1) range(xlims) range(ylims)], 'EdgeColor', 'none', 'FaceColor', faceColor);
                    end
                    
                    plot(axesHandle, time, squeeze(gridData(propIdxs(i), catIdxs(j), :)), 'Color', 'k');
                    
                    set(axesHandle, 'XTick', [], 'YTick', []);
                    set(axesHandle, 'Color', 'none', 'Box', 'off', 'Visible', 'off');
                    xlim(axesHandle, xlims);
                    ylim(axesHandle, ylims);
                    
                    if j == 1
                        ylims = ylim(axesHandle);
                        for k = 1:numel(yticks)
                            annotation(figHandle, 'textbox', [x-0.006 - xMargin*w  (y-i*h + h*(yticks(k)-ylims(1)) / (ylims(end)-ylims(1)))-0.005 0  0], ...
                                'String', ytickLabels{k}, ...
                                'FontSize', fontSizeSmall, ...
                                'EdgeColor', 'none', ...
                                'Margin', 0, ...
                                'HorizontalAlignment', 'Right', ...
                                'VerticalAlignment', 'Middle');
                            annotation(figHandle, 'line', [x  x-0.004] - xMargin*w, (y-i*h + h*(yticks(k)-ylims(1)) / (ylims(end)-ylims(1))) * [1 1] , 'Color', 'k');
                        end
                        annotation(figHandle, 'line', [x; x] - xMargin*w, [
                            y - i*h + yMargin*h + ySpan*h
                            y - i*h + yMargin*h
                            ] , 'Color', 'k');
                    end
                end
            end
            
            %time axis
            axesHandle = subplot('Position', [x+xMargin*w  y-H-yMargin*h  w*xSpan  1e-6]);
            xlim(axesHandle, [0 time(end)]);
            set(axesHandle, 'XTick', [0 10], 'YTick', [], 'TickDir', 'out', 'TickLen', [0.04 0.04], 'FontSize', fontSizeSmall, 'YColor', get(figHandle, 'Color'));
            xlabel(axesHandle, 'Time (h)', 'FontSize', fontSizeLarge);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xlabelPos(1) xlabelPos(2)-50 xlabelPos(3)]);
            tick2text(axesHandle, 'axis', 'x', 'xtickoffset', 7);
            xTicks = getappdata(axesHandle, 'XTickText');
            set(xTicks, 'FontSize', fontSizeSmall, 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            
            %grid
            for i = 0:nProps
                if i > 0
                    txtHandle = annotation(figHandle, 'textbox', [x-0.1-0.008-0.025+0.02  y - (i - 1 / 2) * h .2 0], ...
                        'String', strjoin(' ', propLabels{propIdxs(i), 1}{:}), ...
                        'HorizontalAlignment', 'right', ...
                        'VerticalAlignment', 'Middle', ...
                        'FontSize', fontSizeLarge, ...
                        'EdgeColor', 'none', 'Margin', 0, ...
                        'Interpreter', 'tex');
                    set(txtHandle, 'Position', [x-0.2-0.008-0.018  y - (i - 1 / 2) * h 0.2 0])
                end
                annotation(figHandle, 'line', x + [0 W], y - i * h * [1 1], 'Color', [0.75 0.75 0.75]);
            end
            txtHandles = zeros(nCats, 1);
            for i = 0:nCats
                if i > 0
                    txtHandles(i) = annotation(figHandle, 'textbox', [x+(i-1/2)*w  y+0.03  0  0], ...
                        'String', catLabels{catIdxs(i), 1}, ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'bottom', ...
                        'FontSize', fontSizeMed, ...
                        'EdgeColor', 'none', ...
                        'Margin', 0, ...
                        'FitBoxToText', 'on');
                end
                if i > 0 && catLabels{catIdxs(i), 2} ~= SingleGeneDeletions.WILD_TYPE
                    nGenes = 0;
                    nCorrect = 0;
                    for j = 1:numel(groupCats)
                        if groupCats(j) == catIdxs(i)
                            nGenes = nGenes + ...
                                sum(isGeneImplemented & geneClasses == catLabels{j, 2});
                            nCorrect = nCorrect + ...
                                sum(expEss == modelEss & isGeneImplemented & geneClasses == catLabels{j, 2});
                        end
                    end
                    
                    tmp = strrep(gridMetaData(catIdxs(i)).gene, '_', '');
                    if showGeneSymbols
                        tmp = sprintf('{\\it{%s}}', geneSymbols.(tmp));
                    end
                    
                    txtHandles(i) = annotation(figHandle, 'textbox', [x+(i-1/2)*w  y+0.005  0  0], ...
                        'String', sprintf('(%d, %s)', nGenes, tmp), ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'bottom', ...
                        'FontSize', fontSizeSmall, ...
                        'EdgeColor', 'none', ...
                        'Margin', 0, ...
                        'FitBoxToText', 'on');
                end
                annotation(figHandle, 'line', x + i * w * [1 1], y - [0 H], 'Color', [0.75 0.75 0.75]);
            end
            
            %dendrogram
            annotation(figHandle, 'line', x + ( 1.5) * w * [1 1] + [0 6*w], [y1 y1],  'Color', 'k');
            annotation(figHandle, 'line', x + ( 1.5) * w * [1 1], [y1 y3],  'Color', 'k');
            annotation(figHandle, 'line', x + ( 3.5) * w * [1 1], [y1 y2+0.03],  'Color', 'k');
            annotation(figHandle, 'line', x + ( 6) * w * [1 1], [y1 y2+0.03],  'Color', 'k');
            annotation(figHandle, 'line', x + (7.5) * w * [1 1], [y3 y1],  'Color', 'k');
            
            annotation(figHandle, 'line', x + ( 2.5) * w * [1 1], [y3 y2],  'Color', 'k');
            annotation(figHandle, 'line', x + ( 2.5) * w * [1 1] + [0 2*w], [y2 y2],  'Color', 'k');
            annotation(figHandle, 'line', x + ( 4.5) * w * [1 1], [y2 y3],  'Color', 'k');
            
            annotation(figHandle, 'line', x + (5.5) * w * [1 1], [y3 y2],  'Color', 'k');
            annotation(figHandle, 'line', x + (5.5) * w * [1 1] + [0 w], [y2 y2],  'Color', 'k');
            annotation(figHandle, 'line', x + (6.5) * w * [1 1], [y2 y3],  'Color', 'k');
            
            txtBox = [
                annotation(figHandle, 'textbox', [x+4.5*w-0.15  y1+0.02 0.3 0], 'String', 'Essential')
                annotation(figHandle, 'textbox', [x+3.5*w-0.15  y2+0.02 0.3 0], 'String', 'Macromolecule synthesis')
                annotation(figHandle, 'textbox', [x+6*w-0.15  y2+0.02 0.3 0],   'String', 'Cell cycle')
                ];
            set(txtBox, 'FontSize', fontSizeLarge, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Margin', 0, 'EdgeColor', 'none', 'FitBoxToText', 'on');
            
            %% standarize tick lengths
            paperSize = get(figHandle, 'PaperSize');
            figW = paperSize(1);
            figH = paperSize(2);
            axesHandles = findobj(figHandle, 'type', 'axes', '-and', 'tag', '');
            axesPos = cell2mat(get(axesHandles, 'position'));
            tickLen = 0.0075 * max(max(axesPos(:, 3:4), [], 1) .* [figW figH]);
            for i = 1:numel(axesHandles)
                maxVal = max(axesPos(i, 3:4) .* [figW figH]);
                set(axesHandles(i), 'TickLen', [tickLen / maxVal  0.0250]);
            end
            
            %% save figure
            print(figHandle, outFileName_Figure, '-dpdf', '-rgb');
            close(figHandle);
            
            %% save
            colLabels = {'Class' 'Name' 'No Genes' 'No Simulations' 'Representative Gene Locus Tag' 'Representative Gene Symbol' 'Representative Simulation Batch' 'Representative Simulation Index'};
            content = cell(0, numel(colLabels));
            for i = 1:numel(catIdxs)
                nGenes = 0;
                nSimulations = 0;
                for j = 1:numel(groupCats)
                    if groupCats(j) == catIdxs(i)
                        nGenes = nGenes + gridMetaData(j).nGenes;
                        nSimulations = nSimulations + gridMetaData(j).nSimulations;
                    end
                end
                
                locusTag = gridMetaData(catIdxs(i)).gene;
                if ~isempty(locusTag)
                    symbol = geneSymbols.(strrep(locusTag, '_', ''));
                else
                    symbol = [];
                end
                
                content = [content; {
                    catLabels{catIdxs(i), 1} ...
                    SingleGeneDeletions.DELETION_STRAIN_CLASSES{catLabels{catIdxs(i), 2}, 1} ...
                    nGenes ...
                    nSimulations ...
                    locusTag ...
                    symbol, ...
                    gridMetaData(catIdxs(i)).simGroup ...
                    gridMetaData(catIdxs(i)).simIdx ...
                    }];
            end
            
            if ispc && exist(outFileName_Table, 'file')
                delete(outFileName_Table)
            end
            PrintUtil.printToFile(content, colLabels, outFileName_Table, 'Single Gene Deletions');
        end
        
        
        
    end
    
end

