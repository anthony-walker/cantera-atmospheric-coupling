# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='limonene',
    reactive=True,
    structure=SMILES("CC1=CCC(CC1)C(=C)C"),
)

species(
    label='n-decane',
    reactive=True,
    structure=SMILES("CCCCCCCCCC"),
)

species(
    label='i-octane',
    reactive=True,
    structure=SMILES("CC(C)C(C)C"),
)

species(
    label='toluene',
    reactive=True,
    structure=SMILES("Cc1ccccc1"),
)

species(
    label='benzene',
    reactive=True,
    structure=SMILES("c1ccccc1"),
)

species(
    label='isoprene',
    reactive=True,
    structure=SMILES("C=CC(=C)C"),
)

species(
    label='a-pinene',
    reactive=True,
    structure=SMILES("CC1=CCC2CC1C2(C)C"),
)

species(
    label='b-pinene',
    reactive=True,
    structure=SMILES("CC1(C2CCC(=C)C1C2)C"),
)

species(
    label='farnesane',
    reactive=True,
    structure=SMILES("CCC(C)CCCC(C)CCCC(C)C"),
)

species(
    label='H2S',
    reactive=True,
    structure=SMILES("S"),
)
species(
    label='O2',
    structure=SMILES("[O][O]"),
)

species(
    label='N2',
    reactive=False,
    structure=SMILES("N#N")
)

species(
    label='AR',
    reactive=False,
    structure=SMILES("[Ar]")
)

# Reaction systems
simpleReactor(
    temperature=(2000,'K'),
    pressure=(5.0,'bar'),
    initialMoleFractions={ "n-decane":5.195e-03, "i-octane":4.020e-03, "toluene":2.960e-03, "H2S":1.522e-07, "farnesane":7.610e-04, "limonene":7.610e-04, "a-pinene":7.610e-04, "b-pinene":7.610e-04, "N2":7.779e-01, "O2":2.069e-01
    },
    terminationTime=(1e1,'s'),
)

simpleReactor(
    temperature=(1500,'K'),
    pressure=(5.0,'bar'),
    initialMoleFractions={ "n-decane":5.195e-03, "i-octane":4.020e-03, "toluene":2.960e-03, "H2S":1.522e-07, "farnesane":7.610e-04, "limonene":7.610e-04, "a-pinene":7.610e-04, "b-pinene":7.610e-04, "N2":7.779e-01, "O2":2.069e-01
    },
    terminationTime=(1e1,'s'),
)

simpleReactor(
    temperature=(1000,'K'),
    pressure=(5.0,'bar'),
    initialMoleFractions={ "n-decane":5.195e-03, "i-octane":4.020e-03, "toluene":2.960e-03, "H2S":1.522e-07, "farnesane":7.610e-04, "limonene":7.610e-04, "a-pinene":7.610e-04, "b-pinene":7.610e-04, "N2":7.779e-01, "O2":2.069e-01
    },
    terminationTime=(1e1,'s'),
)

simpleReactor(
    temperature=(500,'K'),
    pressure=(5.0,'bar'),
    initialMoleFractions={ "n-decane":5.195e-03, "i-octane":4.020e-03, "toluene":2.960e-03, "H2S":1.522e-07, "farnesane":7.610e-04, "limonene":7.610e-04, "a-pinene":7.610e-04, "b-pinene":7.610e-04, "N2":7.779e-01, "O2":2.069e-01
    },
    terminationTime=(1e1,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceMoveToCore=0.01,
    toleranceInterruptSimulation=0.01,
)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2500,'K',8),
    pressures=(0.01,10,'bar', 5),
    interpolation=('Chebyshev', 6, 4),
)

quantumMechanics(
        software='mopac',
        method='pm3',
        fileStore='QMfiles',
        scratchDirectory = None,
        onlyCyclics = True,
        maxRadicalNumber = 0,
        )

options(
    units='si',
    generatePlots=False,
    saveEdgeSpecies=True,
    verboseComments=True
)
