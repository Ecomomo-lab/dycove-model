[SedimentFileInformation]
    FileCreatedBy         = Deltares, FM-Suite DFlowFM Model Version 4.3.0.0, DFlow FM Version 1.2.110.67911M
    FileCreationDate      = Sat Dec 06 2025, 15:43:26
    FileVersion           = 02.00
[SedimentOverall]
    Cref                  = 1600                   [kg/m3]   Reference density for hindered settling calculations
[Sediment]
    Name                  = #sedimentsand#                   Name of sediment fraction
    SedTyp                = sand                             Must be "sand", "mud" or "bedload"
    IniSedThick           = 5                      [m]       Initial sediment layer thickness at bed
    FacDss                = 1                      [-]       Factor for suspended sediment diameter
    RhoSol                = 2650                   [kg/m3]   Specific density
    TraFrm                = -2                               Integer selecting the transport formula
    CDryB                 = 1600                   [kg/m3]   Dry bed density
    SedDia                = 0.0002                 [m]       Median sediment diameter (D50)
    IopSus                = 0                                Option for determining suspended sediment diameter
    Pangle                = 0                      [degrees] Phase lead angle
    Fpco                  = 1                      [-]       Coefficient for phase lag effects
    Subiw                 = 51                               Wave period subdivision
    EpsPar                = False                            Use Van Rijn's parabolic mixing coefficient
    GamTcr                = 1.5                    [-]       Coefficient for grain size effect
    SalMax                = 0                      [ppt]     Salinity for saline settling velocity

[Sediment]
    Name                  = #sedimentorganic#                Name of sediment fraction
    SedTyp                = mud                              Must be "sand", "mud" or "bedload"
    IniSedThick           = 0                      [m]       Organic sediment starts absent; DYCOVE supplies OM
    FacDss                = 1                      [-]       Factor for suspended sediment diameter
    RhoSol                = 1000                   [kg/m3]   Specific density
    TraFrm                = -3                               Cohesive sediment transport formulation
    CDryB                 = 100                    [kg/m3]   Dry bed density used for mass-volume conversion
    SalMax                = 31                     [ppt]     Salinity for saline settling velocity
    WS0                   = 0.00025                [m/s]     Settling velocity in fresh water
    WSM                   = 0.00025                [m/s]     Settling velocity in saline water
    EroPar                = 0.0001                 [kg/m2/s] Erosion parameter
    TcrSed                = 1000                   [N/m2]    Critical stress for sedimentation
    TcrEro                = 0.2                    [N/m2]    Critical stress for erosion
    TcrFluff              = 4.94065645841247E-324  [N/m2]    Critical stress for fluff-layer erosion
