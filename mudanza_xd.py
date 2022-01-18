# -*- coding: utf-8 -
# Mudanza_xD.ipynb


# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >


''' Paquetes necesarios '''

from unittest import signals
import  numpy               as  np
import  matplotlib.pyplot   as  plt
import  bioread
from    scipy.signal import savgol_filter, find_peaks, find_peaks_cwt
# agregar find peaks


# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >


def open_acq( acq_name ):
    try:
        data = bioread.read( acq_name )
        return data
    except:
        print ( "\n\t [!] Error al abrir el el documento. \n" )


# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >


'''
busca los indices donde se encuentra el canal de referencia -> tupla
'''

def get_indexes( reference ):
    return np.nonzero( reference == 5 )


# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >


''' 
Esta funcion es particularmente importante. Lo que hace es guardar las señales de
respuesta con el largo especificado. Una es un arreglo de una dimensión con las
señales de corrido.  El otro es una lista, que guarda en cada elemento el arreglo
de una sola señal.
'''

def get_response ( register, reference, RESPONSE_SIZE ):

    response = np.array( [] )   # arreglo de tipo numpy
    l_response = []             

    # i_reference es una tupla y contiene un numpy.ndarray
    i_reference = get_indexes( reference )

    i = 0; j = 0
    for index in i_reference[0]:
        if j > 0 and i_reference[0][j] == i_reference[0][j-1] + 1:
            j += 1
            continue
        else:
            # registro
            # if register[index : index + RESPONSE_SIZE].max() > 0.5:
            i += 1; j += 1
            l_response.append( register[index : index + RESPONSE_SIZE ] )
            response = np.append( response, register[ index : index + RESPONSE_SIZE ] )  

    return l_response, response


# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >


def filter( l_response, response ):

    response_fitted = np.array( [] )
    l_response_fitted = [] 

    window = 33
    grade = 3
    # se aplica el filtro sólo a las respuestas
    response_fitted = savgol_filter( response, window, grade, mode='nearest')
    for elem in l_response:
        l_response_fitted.append(savgol_filter( elem, window, grade, mode='nearest'))
    
    return l_response_fitted, response_fitted


# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >


def write_response ( l_response, nombre ):

    encabezado = ''
    for elem in range( len( l_response ) ):
        if elem > 0 and elem < len( l_response ):
            encabezado += '|'
        encabezado += "Respuesta" + ' '
        encabezado +=  str( elem + 1 ) 

    try:
        np.savetxt( nombre + '.txt', np.array( l_response ).T, delimiter = "|", header = encabezado, comments = '', fmt='%.15f' )
        print( '\n\t [*] El archivo ' + nombre + '.txt' + ' se generó exitosamente. \n' )
    except:
        print ( '\n\t [!] Error al crear ' + nombre + '.txt \n' )
        

# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >


def diff_mean( l_response_fitted ):

    response_diff_m = []
    l_response_diff_m = []

    for response in l_response_fitted:
        response_diff_m = []
        for i in range( len( response ) ):
            if i > 0 and i < len( response ) - 1:
                response_diff_m.append( ( response[i+1] - response[i-1]) / 2 )
            else:
                response_diff_m.append( 0 )

        l_response_diff_m.append(response_diff_m)

    return l_response_diff_m


# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >


def invertir( lista ):
    
    for i in range(len(lista)):
        lista[i] = -lista[i]
        
    return lista


# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >


def get_peaks( l_response_fitted ):

    l_index = []

    for signal in l_response_fitted:
        i_positive, _ = find_peaks(signal, prominence = 0.5)
        inv_response = np.copy(signal)
        i_negative, e = find_peaks(invertir(inv_response), prominence = 0.5)
        index = np.concatenate((i_positive, i_negative))
        index = np.sort(index)
        l_index.append( index )
    
    return l_index 


# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >

'''
Esta función recorré los índices de cada señal, pues, estos son calculados en cada señal
pero en las respuestas concatenadas debe irse separando cada RESPONSE_SIZE valores.
'''

def stretch_index( l_peaks, RESPONSE_SIZE ):
    n = 0
    l_index = np.array([])
    for arreglo in l_peaks:
        for elem in arreglo:
            l_index = np.append( l_index, elem + n * RESPONSE_SIZE)
        n += 1

    return l_index.astype(int)


# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >

def fill_peaks ( l_peaks_alone ):

    for i in range( len(l_peaks_alone) ):
        length = len (l_peaks_alone[i])
        if length == 0:
            l_peaks_alone[i] = np.concatenate( (l_peaks_alone[i], np.array([0,0,0,0])) )
        elif length == 1:
            l_peaks_alone[i] = np.concatenate( (l_peaks_alone[i], np.array([0,0,0])) )
        elif length == 2:
            l_peaks_alone[i] = np.concatenate( (l_peaks_alone[i], np.array([0,0])) )
        elif length == 3:
            l_peaks_alone[i] = np.concatenate( (l_peaks_alone[i], np.array([0])) )

    for elem in l_peaks_alone:
        print (elem)
    
    return l_peaks_alone


# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >


def plot_response_IO_artefact( response, response_fitted, l_peaks, RESPONSE_SIZE, acq_name ):

    l_index = stretch_index( l_peaks, RESPONSE_SIZE )

    n = 0
    for elem in l_peaks:
        print(response[elem[0] + n * RESPONSE_SIZE], response[elem[1] + n * RESPONSE_SIZE ], response[elem[2]+ n * RESPONSE_SIZE], response[elem[3]+ n * RESPONSE_SIZE])
        n+=1


    plt.figure( acq_name, facecolor="#F1F2F3" )
    plt.grid( axis = 'y', color = 'gray', linestyle = "dotted", linewidth = 0.4 )
    plt.grid( axis = 'x', color = 'gray', linestyle = "dotted", linewidth = 0.4 )
    plt.box()

    plt.plot( response, linewidth = 1 )
    plt.plot(response_fitted, linewidth = 0.7, color = "red" )
    plt.scatter( l_index, response[l_index], color = "orange", marker = "x")

    print ( "\n\t [*] El archivo grafica.jpg se generó exitosamente. \n" )

    plt.show()


# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >


''' main de pruebas '''

def main ():

    RESPONSE_SIZE = 5000

    acq_name = "/home/luuisraam/Escritorio/Mudanza/R54 Ep + Lev SE 21 noviembre 2018 Registro 8 abril 2019 HD GD Señal 5 Curva IO Pre Lev.acq"

    # -- [1] Apertura del archivo acq y lectura de los datos
    data = open_acq( acq_name )
    register = data.channels[0].data    # Es el canal de respuesta
    reference = data.channels[1].data   # Es el canal de referencia

    # -- [2] cortado y filtrado de las señales
    l_response, response = get_response ( register, reference, RESPONSE_SIZE )
    l_response_fitted, response_fitted = filter( l_response, response )

    print ("\n\t  Se encontraron " + str(len(l_response)) + " respuestas. \n")  
    
    # -- [3] encuentra los maximos y minimos

    l_peaks_alone = get_peaks ( l_response_fitted )
    l_peaks = fill_peaks ( l_peaks_alone )
    
    # Busca las derivadas

    l_response_diff_m = diff_mean( l_response_fitted )

    # -- [4] Genera los archivos

    write_response( l_response, acq_name + "_Respuesta" )
    write_response( l_response_fitted, acq_name + "_Respuesta_suave" )


    plot_response_IO_artefact( response, response_fitted, l_peaks, RESPONSE_SIZE, acq_name )

    return 0


main ()
