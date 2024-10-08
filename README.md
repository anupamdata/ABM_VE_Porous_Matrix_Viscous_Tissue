# ABM_VE_Porous_Matrix_Viscous_Tissue
This code simulates the growth of viscous tissue encapsulated in a porous viscoelastic matrix.

Agent-Based Model for Viscous Tissue Growth in a Porous Matrix:

This code, developed by Ms. Kirti Kashyap and Dr. Anupam Gupta, simulates the growth of viscous tissue encapsulated in a porous viscoelastic matrix using an agent-based model. It builds upon the work presented in https://github.com/anupamdata/ABM_VE_Matrix_Viscous_Tissue.git (further details: https://doi.org/10.1038/s41563-022-01400-4) by incorporating the porosity of the extracellular matrix (ECM).

Features:

Modular Fortran Code: Written in Fortran for compatibility with various Fortran compilers.

User-Friendly Input: The para.in file allows you to specify parameters for simulating tissue growth.

Porosity Control: Simulate the impact of ECM porosity on tissue growth.

Key Parameters:

mu (line 19): Controls the viscosity of the matrix.

Es (line 15): Controls the elasticity of the matrix.

gm (line 18): Controls the viscosity of the tissue.

poro_phi (line 9): Controls the fraction of passive beads (Degree of porosity).

Sample Usage:

The provided sample parameters in para.in are configured to simulate a branching tissue growing within a 30% porous matrix.
