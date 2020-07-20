# continentalDeformationMAT
A MATLAB package for modelling geologic scale continental deformation, using the thin viscous sheet approximation 

# Motivation
Continental collision produces some of the most striking tectonic features on the Earth's surface. For example, collision of the Indian tectonic plate into Asia at 60 Ma caused the formation of the Himalayan range and surrounding Tibetan Plateau, the highest elevation on Earth. These tectonic processes caused huge amounts of uplift and crustal thickening over a far horizontal extent. To better understand these processes, numerical models can be used to test mechanisms and physical properties. 

Here, I present MATLAB package for modelling geologic scale continental deformation, using the thin viscous sheet approximation of [England and McKenzie (1982)](https://doi.org/10.1111/j.1365-246X.1982.tb04969.x). Appropriate corrections and changes are made according to [later publications](https://doi.org/10.1111/j.1365-246X.1983.tb03328.x). 

# Theory
The idea behind the thin viscous sheet approximation is to represent the crust as a thin layer that is isostatically in equilibrium with the lithosphere, and that velocities within the crust do not vary with depth. This flow is incompressible and is driven by tractions and the gradients in crustal thickness. Vertical gradients of deviatoric stresses are also neglected. Considering vertical gradients in forces and rheology, it can be shown that pressure in a column of crust is

<img src="https://render.githubusercontent.com/render/math?math=p=\tau_{zz}-\int_{0}^{z}\rho g\ dz'.">

Consider a level at the base of the lithosphere for isostatic balance (z=0). The vertical average of this pressure through the lithosphere and any topography gives

<img src="https://render.githubusercontent.com/render/math?math=\bar{p}=\bar{\tau}_{zz}-\frac{g \rho_c}{2L}(1-\rho_c/\rho_m)S^2+\frac{g \rho_m L}{2}.">

The term <img src="https://render.githubusercontent.com/render/math?math=\bar{\tau}_{zz}"> is calculated through the constitutive equation

<img src="https://render.githubusercontent.com/render/math?math=\bar{\tau}_{ij}=B\dot{E}^{(1/n-1)}\dot{\varepsilon}_{ij}.">

This describes a power-law rheology, where n=1 is a Newtonian fluid. The parameter B is a constant which encompasses depth averaged viscosity, and <img src="https://render.githubusercontent.com/render/math?math=\dot{E}"> is the second invariant of the strain rate tensor.

Substituting <img src="https://render.githubusercontent.com/render/math?math=\bar{p}"> and <img src="https://render.githubusercontent.com/render/math?math=\bar{\tau}_{zz}"> into force-balance equations in the horizontal directions gives

<img src="https://render.githubusercontent.com/render/math?math=\nabla^2 u=2Ar\dot{E}^{(1-1/n)}S\nabla S-3\nabla(\nabla\cdot u)-2(1/n-1)\dot{E}^{-1}\nabla\dot{E}\cdot\dot{\varepsilon}-2(1/n-1)\dot{E}^{-1}(\nabla\cdot u)\nabla\dot{E}.">

In this representation, only horizontal derivatives are considered. The terms are non-dimensional and the Argand number (Ar) is defined as:

<img src="https://render.githubusercontent.com/render/math?math=Ar = \frac{g\rho_c(1-\rho_c/\rho_m)L}{B(u_0/L)^{1/n}}.">

This number describes the relative contributions of forces arising from variations in crustal thicknesses and the force required to deform a fluid. For large Ar, the crustal gravitational forces succumb to deformation - the lithosphere is weak (a cheese analogue might be a ripe brie). Small Ar flows resist deformation and changes in crustal thickness, forming rigid blocks (think feta).

To constrain time-dependence, the continuity equation is written in the form:

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial S}{\partial t}=- \nabla\cdot (Su).">

# Implementation and computation

Equations are solved on a NxN grid (e.g. 32x32), through numerical approximations to the equations for velocity and crustal thickness. The derivative terms in the RHS of the velocity equation are approximated through centre difference schemes. To solve for the velocity field on the LHS, the interior nodes of the RHS are initially set to zero and the new velocity field is inverted for using a Poisson equation solving routine.

<img src="https://user-images.githubusercontent.com/38541020/87995475-74b9c380-caa4-11ea-9c47-e6c0b5bd77ac.png" width="500" height="auto"/>

Boundary conditions are applied to the velocity fields. These new velocity fields are then used in the centre difference schemes to calculate a new RHS, which is solved again to attain a new solution for the velocities (again applying boundary conditions). This new solution is applied to the old solution until the solution converges. For calculations, very small convergence parameters must be used.

To advance the solution in time, as well as calculate changes in crustal thickness, numerical approximations to the continuity equation are used. The LHS is determined using forward differences, and the RHS is found using an 'upwind' scheme, as the form of the equation lends similarities to a material derivative. The timestep is chosen to satisfy the Courant-Friedrichs-Lewy criterion.

<img src="https://user-images.githubusercontent.com/38541020/87989771-7d57cd00-ca97-11ea-967f-e772b12c35c3.png" width="500" height="auto"/>

# Examples

To illustrate the uses of this package I show the results of two simulations: one at Ar=1; and another at Ar=10. Both simulations are Newtonian (n=1) and are solved over a range of 5 million years. Diagnostic parameters are shown below for each model ()

<img src="https://user-images.githubusercontent.com/38541020/87987748-0240e780-ca94-11ea-9041-29bdd6645537.png" width="500" height="auto"/><img src="https://user-images.githubusercontent.com/38541020/87987782-0cfb7c80-ca94-11ea-87b6-7965e35b49d0.png" width="500" height="auto"/>

<img src="https://user-images.githubusercontent.com/38541020/87987422-7c24a100-ca93-11ea-9592-a0246925571b.png" width="500" height="auto"/><img src="https://user-images.githubusercontent.com/38541020/87987473-8d6dad80-ca93-11ea-9b2d-4c75dbcff4ab.png" width="500" height="auto"/>

# Examples

<img src="https://user-images.githubusercontent.com/38541020/87986639-3f0bdf00-ca92-11ea-9e81-d23afbbff34f.gif" width="500" height="auto"/><img src="https://user-images.githubusercontent.com/38541020/87987207-19cba080-ca93-11ea-9995-9f8749467be7.gif" width="500" height="auto"/>

