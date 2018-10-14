var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#MRIReco.jl-1",
    "page": "Home",
    "title": "MRIReco.jl",
    "category": "section",
    "text": "Magnetic Resonance Imaging Reconstruction"
},

{
    "location": "index.html#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "MPIReco is a Julia packet for magnetic resonance imaging. It contains algorithms for the simulation and reconstruction of MRT data and is both easy to use and flexibly expandable.Both direct and iterative methods are available for image reconstruction. In particular, modern compressed sensing algorithms such as ADMM can be used.The MRT imaging operator can be set up for a variety of scanning patterns (cartesian, spiral, radial, ...) and can take into account field inhomogeneity as well as the use of coil arrays. The operator can be quickly evaluated using NFFT-based methods.One strength of the package is that it is strongly modular and uses high quality Julia packages. These are e.g.NFFT.jl and FFTW.jl for fast Fourier transformations\nWavelets.jl for sparsification\nLinearOperators.jl in order to be able to divide the imaging operator modularly into individual parts\nLinearSolver.jl for modern algorithms for solving linear optimization problemsThis interaction allows new algorithms to be easily integrated into the software framework. It is not necessary to program in C/C++ but the advantages of the scientific high-level language Julia can be used."
},

{
    "location": "index.html#Status-1",
    "page": "Home",
    "title": "Status",
    "category": "section",
    "text": "MRIReco.jl is work in progress and in some parts not entirely optimized. In particular the FFT and NFFT implementation are currently limited to the CPU and do not support GPU acceleration yet."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Start julia and open the package mode by entering ]. The enteradd https://github.com/tknopp/LinearSolver.jl\nadd https://github.com/MagneticResonanceImaging/MRIReco.jlThis will install the packages LinearSolver.jl, MRIReco.jl, and all its dependencies."
},

{
    "location": "gettingStarted.html#",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "page",
    "text": ""
},

{
    "location": "gettingStarted.html#Getting-Started-1",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "section",
    "text": "We will start with a very simple example and perform simple simulation and reconstruction based on a shepp logan phantom. The programm looks like this# image\nN = 256\nx = shepp_logan(N)\n\n# simulation\nparams = Dict{Symbol, Any}()\nparams[:simulation] = \"nfft\"\nparams[:trajName] = \"Cartesian\"\nparams[:numProfiles] = floor(Int64, N)\nparams[:numSamplingPerProfile] = N\n\naqData = simulation(x, params)\naqData = MRIReco.sample_kspace(aqData, redFac, \"poisson\", calsize=5)\naqData = convertUndersampledData(aqData)\n\n# reco\nparams[:reco] = \"simple\"    # encoding model\nparams[:shape] = (N,N)\nparams[:sparseTrafoName] = \"Wavelet\" #sparse trafo\nparams[:regularization] = \"L1\"       # regularization\nparams[:lambdL1] = 1.e-3\nparams[:solver] = \"admm\"    # solver\nparams[:iterations] = 1000\nparams[:ρ] = 1.0e-1\n\nx_reco = reconstruction(aqData, params)We will go through the program step by step. First we create a 2D shepp logan phantom of size N=256. Then we setup a dictionary that defines the simulation parameters. Here, we chose a simple Cartesian trajectory ..."
},

]}
