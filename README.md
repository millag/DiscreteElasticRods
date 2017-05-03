This project is a stand alone hair simulation tool written as part of my master thesis at BU (link here) and incorporates a few components - 
1 acurately capturing the dynamics of a single hair strand of an arbitrary shape - straight, wavy or curled
2 coarse approximation of hair to hair interactions
3 OpenGL preview of the simulated data
4 export of the hair system as standard .obj files to import in other DCC software as f.e. Houdini(see below).

Single strand dynamics was implemented with the intention to capture all the elastic energy contributions - streching, bending and twisting,
where stretching and bending are modeled by the deformations of the space curve, while twisting - by the rotation of a material frame associated 
with each point on the curve. Implementation is based on the paper "Discrete Elastic Rods" by Bergou et al.,
ACM Transactions of Graphics (SIGGRAPH) 2008 (link here), with a few modifications made to the way coupling rods to rigid body and enforcing 
rod inextensibility is handled. Position-based dynamics is used instead of the proposed manifold projection method. This achieves better performance
and helps improve the stability of the simulation.

OpenGL preview:
https://vimeo.com/104214025

Imported in Houdini and rendered with Mantra:
Video:
https://vimeo.com/107726073
Images:
![curly1](https://cloud.githubusercontent.com/assets/133900/25644666/8deb69e0-2fb0-11e7-8a5e-b12966618a28.jpg)
![curly2](https://cloud.githubusercontent.com/assets/133900/25644672/9255b5bc-2fb0-11e7-89de-55003e8e432a.jpg)
![wavy1](https://cloud.githubusercontent.com/assets/133900/25644673/961af5c2-2fb0-11e7-88ac-8ae2eac4a33b.jpg)
