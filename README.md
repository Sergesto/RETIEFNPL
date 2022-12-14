## RETIEFNPL
### RETIculado Elasticidad Finita No lineal PLástica
#### para MATLAB o GNU–Octave

- Estructuras Reticuladas 3D
- Elasticidad finita, material de Saint Venant–Kirchhoff
- Deformación de Green–Lagrange, hipótesis de pequeñas deformaciones
- Plasticidad, algoritmo de mapeo de retorno

### Contenido
- retiefnpl.m (main program)
- amr.m / Algoritmo de mapeo de retorno –function–
- matrixfn.m / Matriz Kr de la estructura (matriz reducida por condiciones de borde) –function–
- gradconj.m / Método del gradiente conjugado –function–

### Referencias
- Simo, J.C. and Hughes, T.J.R. (1998) Computational Inelasticity. Springer-Verlag, New York.
- Bazzano, J.B., Pérez Zerpa, J. (2017) Introducción al Análisis No Lineal de Estructuras, Facultad de Ingeniería, Universidad de la República, Montevideo, Uruguay.
- Morton E. Gurtin (1981) Introduction to Continuum Mechanics (Mathematics in Science and Engineering, Volume 158), Carnegie-Mellon University, Pittsburgh, Pennsylvania.
- A. Canelas (2013) Elasticidad Finita –apuntes del curso–, Facultad de Ingeniería, Universidad de la República, Montevideo, Uruguay.

### Ejemplo

![fig1](https://user-images.githubusercontent.com/104937664/182744989-27ebfd10-c7e4-4781-80af-2d0f59a2e547.jpg)
![fig2](https://user-images.githubusercontent.com/104937664/182745127-7b927523-817e-400a-abb3-d0db8ce2ea0c.jpg)
![fig3](https://user-images.githubusercontent.com/104937664/182745158-4e29bd12-42a3-4f43-a361-6544ef470222.jpg)
