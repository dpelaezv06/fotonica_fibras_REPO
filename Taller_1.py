import numpy as np
import modos_opticaRayos as ray

''' parametros del guia de onda '''
n_core = 1.5 #indice de refraccion del nucleo
n_cleavy = 1 #indice de refraccion del recubrimiento
n_substract = 1 #indice de refraccion del sustrato
espesor = 1 #espesor del guia de onda
longitud_onda = 1 #longitud de onda usada en la fibra optica

''' calculo de los valores solucion de la ecuacion trascendente'''
for modo in range(0, 3): #se recorren los valores de modos desde el 0 hasta el 3
    angulo_optimo, valor_min = ray.optimizar_TERayos(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda) #calculo del angulo para el modo TM
    ''' printeamos por consola el valor del angulo y de la ecuacion trascendednte en ese angulo. Si el valor de la ecuacion no es cero, entonces ese modo no es 
    permitido en el guia de onda '''
    indice_efectivo = ray.n_efectivo(n_core, angulo_optimo)
    print(f"Modo TE {modo}: Angulo optimo = {angulo_optimo:.6f}, n efectivo = {indice_efectivo:.6f}") #se arroja el print de los datos


    angulo_optimo, valor_min = ray.optimizar_TMRayos(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda) #calculo del angulo para el modo TE
    indice_efectivo = ray.n_efectivo(n_core, angulo_optimo)
    ''' printeamos por consola el valor del angulo y de la ecuacion trascendednte en ese angulo. Si el valor de la ecuacion no es cero, entonces ese modo no es 
    permitido en el guia de onda '''
    print(f"Modo TM {modo}: Angulo optimo = {angulo_optimo:.6f}, n efectivo = {indice_efectivo:.6f}") #se arroja el print de los datos

