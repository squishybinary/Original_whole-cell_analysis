function [simulation, states, options, parameters, fittedConstants, metadata] = processRawData(directoryOrWID, fileName, timeInit, timeFin, downsampleStepSec, downsampleType)
    %Initialize
    setWarnings();
    setPath();
    setPreferences();

    %Import classes
    import edu.stanford.covert.cell.sim.Simulation;
    import edu.stanford.covert.cell.sim.util.DatabaseLogger;
    import edu.stanford.covert.cell.sim.util.DiskLogger;
    import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
    import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
    import edu.stanford.covert.db.MySQLDatabase;

    %Load metadata
    if ischar(directoryOrWID)
        directoryOrWID = SimulationDiskUtil.getSimulation(directoryOrWID);
        metadata = DiskLogger.loadMetadata(directoryOrWID);
    else
        if ~exist('database', 'var')
            database = MySQLDatabase(getConfig());
        end
        metadata = DatabaseLogger.loadMetadata(database, directoryOrWID);
    end

    %Construct simulation object
    simulation = CachedSimulationObjectUtil.load();
    %metadata.revision > [], metadata or ()

    %Load options, parameters, fitted constants based on data stored in disk/database
    if ischar(directoryOrWID)
        if nargin < 3
            timeInit = [];
        end
        if nargin < 4
            timeFin = [];
        end
        if nargin < 5
            downsampleStepSec = [];
        end
        if nargin < 6
            downsampleType = [];
        end
        [states, metadata, options, parameters, fittedConstants, randStreamStates] = ...
            DiskLogger.load(directoryOrWID, '-independent', timeInit, timeFin, downsampleStepSec, downsampleType);
    else
        if ~exist('database', 'var')
            database = MySQLDatabase(getConfig());
        end
        if nargin > 1
            warning('WholeCell:warning', 'Additional options ignored');
        end
        [states, metadata, options, parameters, fittedConstants, randStreamStates] = ...
            DatabaseLogger.load(simulation, database, directoryOrWID);
    end

    %Apply options, parameters, constants, time courses to simulation
    simulation.applyOptions(options);
    simulation.applyParameters(parameters);
    simulation.applyFittedConstants(fittedConstants);
    simulation.applyRandStreamStates(randStreamStates);
    simulation.allocateMemoryForState(numel(states.Time.values));

    for i = 1:numel(simulation.states)
        state = simulation.states{i};
        stateID = state.wholeCellModelID(7:end);
        for j = 1:numel(state.stateNames)
            name = state.stateNames{j};
            try
                state.(name) = states.(stateID).(name);
            catch %#ok<CTCH>
                warning('WholeCell:warning', 'Data not provided for %s %s state', state.name, name);
            end
        end
    end

    for i = 1:numel(simulation.processes)
        process = simulation.processes{i};
        process.copyFromState();
    end

    simulation.applyPerturbationsToConstants();

    %Store simulation object in data folder / cache
    %pass inoutdir / filename
    CachedSimulationObjectUtil.store(simulation, [], fileName);
    % (simulation, [] OR kbWid, fileName or this.simCache)

    %Cleanup
    if exist('database', 'var')
        database.close();
    end
end

