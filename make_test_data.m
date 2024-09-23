function testdata = make_test_data(scenario,n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% location of random sample
testdata.lat= repmat(30,n,1);
testdata.lon= repmat(10,n,1);
testdata.altitude=repmat(500,n,1);

switch scenario
    case 'step'
        testdata.t = [1500];                         % step change timing
        testdata.e = [20,50,100,300,50,400,50,...    % old erosion rates of different catchments
                      1000,2000,4000,200,50,100,500];% new erosion rates of catchments
        testdata.changeVariable = [];                % this is variable for convinience, so that I dont need a swith between loss/change or no change factor later on
    case 'samestep'
        testdata.t = [1500];                        % step change timing
        testdata.e = [20,50,100,40,60,80,90];       % old erosion rates of different catchments
        testdata.chg = [20];                        % change factor of erosion rate at time t
    case 'samebackground_step'
        testdata.t = [1500];                         % step change timing
        testdata.e = 50;                             % old erosion rates of different catchments
        testdata.chg = [10,20,5,40,2,25,30];       % change factor of erosion rate at time t
    case 'samebackground_samestep'
        testdata.t = 1000;
        testdata.e = 50;
        testdata.chg= 10;
    case 'spike'
        testdata.t =  [1500];                           % soil loss timing
        testdata.e =  [20,50,100,300,30,400,50];        % background erosion rates of different catchments
        testdata.loss =[50,20,1,  30, 10, 20, 40];       % soil loss in cm
    case 'samespike'
        testdata.t = [1500];                            % soil loss timing
        testdata.e = [20,50,100,300,50,400,50];         % background erosion rates of different catchments
        testdata.loss = [50];                            % soil loss in cm
    case 'samebackground_spike'
        testdata.t =   [1500];                          % soil loss timing
        testdata.e =   [50];                            % background erosion rates of different catchments
        testdata.loss = [50,10,100,30,20,50,5];         % soil loss in cm
    case 'samebackground_samespike'
        testdata.t =   [1500,500];                      % soil loss timing
        testdata.e =   [50];                            % background erosion rates of different catchments
        testdata.loss = [10,20];                        % soil loss in cm
end

testdata.steps = length(testdata.t);
if isfield(testdata,'chg')
    testdata.changeVariable = testdata.chg;
elseif isfield(testdata,'loss')
    testdata.changeVariable = testdata.loss;
else
    testdata.changeVariable = [];
end


end