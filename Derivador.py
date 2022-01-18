import  numpy               as  np
import  matplotlib.pyplot   as  plt
import  bioread
from    tkinter import Tk, messagebox
from    tkinter import filedialog as FileDialog
from    scipy.signal import savgol_filter


'''
Funciones:
    search_file: ( None ) -> str
    open_acq: ( None ) -> Datafile | None
    get_indexes: ( numpy.ndarray ) -> Tupla
    get_response ( numpy.ndarray, numpy.ndarray, int ) -> numpy.ndarray, list[numpy.ndarray]
    plot_response: ( list, numpy.ndarray, int, str ) -> None
    write_response: (list) -> None
    main: ()
'''

# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
def search_file( ):
    return FileDialog.askopenfilename(title="Abrir archivo de registros", filetypes = (("ACQ Knowledge","*.acq"),) )


# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
def open_acq( acq_name ):
    try:
        data = bioread.read( acq_name )
        return data
    except:
        print ( "\n\t [!] Error al abrir el el documento. \n" )


# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
def get_indexes( reference ):
    # busca los indices donde se encuentra el canal de referencia -> tupla
    return np.nonzero( reference == 5 )



# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
# Esta funcion es particularmente importante. Lo que hace es guardar las señales de
# respuesta con el largo especificado. Las guarda para plotearla. Una es un arreglo
# de una dimensión con las señales decorrido.  El otro es una lista,  que guarda en
# cada elemento el arreglo de una sola señal.
def get_response ( register, reference, response_size ):

    response = np.array( [] )
    l_response = []

    # i_reference es una tupla y contiene un numpy.ndarray
    i_reference = get_indexes( reference )

    i = 0; j = 0
    for index in i_reference[0]:
        if j > 0 and i_reference[0][j] == i_reference[0][j-1] + 1:
            j += 1
            continue
        else:
            # if register[index : index + response_size].max() > -5:
            i += 1; j += 1
            l_response.append( register[index: index + response_size ] )
            response = np.append( response, register[index:index+response_size] )

    return response, l_response



# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
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



# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
def plot_response_2( l_response, response, response_fitted, response_size, acq_name ):

    plt.figure( "Respuestas encontradas" +  acq_name, facecolor="#F1F2F3" )

    tiempo = []
    t = 0
    for elem in response:
        t = t+1
        tiempo.append(t)
    
    plt.subplot( 2, 1, 1 )
    plt.plot( l_response[1:], linewidth = 0.8 )
    plt.grid( axis = "x", color = "gray", linestyle = "dotted", linewidth = 0.4 )
    plt.grid( axis = "y", color = "gray", linestyle = "dotted", linewidth = 0.4 )
    plt.ylabel ( "mV" )

    plt.subplot( 2, 1, 2 )
    plt.plot( tiempo, response_fitted, linewidth = 0.5, color = "orange" )
    plt.scatter( tiempo, response, s=0.8 )
    plt.box()
    plt.ylabel ( "mV" ); plt.xlabel( "tiempo" )
    plt.grid( axis = 'y', color = 'gray', linestyle = "dotted", linewidth = 0.4 )
    for i in range( 1, len( l_response[0]) ):
        if i % 5 == 0:
            plt.axvline( x = i * response_size, linestyle = "dotted", color = "r", linewidth = 0.8 )
        else:
            plt.axvline( x = i * response_size, linestyle = "dotted", color = "gray", linewidth = 0.5 )
    '''
    
    plt.subplot( 3, 1, 3 )
    plt.scatter( tiempo, response_fitted, s=0.8 )
    plt.box()
    plt.ylabel ( "mV" ); plt.xlabel( "tiempo" )
    plt.grid( axis = 'y', color = 'gray', linestyle = "dotted", linewidth = 0.4 )
    for i in range( 1, len( l_response[0]) ):
        plt.axvline( x = i * response_size, linestyle = "dotted", color = "r", linewidth = 0.8 )
    '''

    # plt.savefig("Grafica.png")
    # print ( "\n\t [*] El archivo grafica.png se generó exitosamente. \n" )

    plt.show()



# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
def diff_right( l_response_fitted ):

    l_response_diff_r = []

    for response in l_response_fitted:
        l_response_diff_r.append( np.diff(response) )
    
    return l_response_diff_r


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



# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
def int_trapezium( l_response_fitted ):

    response_int = []
    l_response_int = []

    for response in l_response_fitted:
        response_diff_m = []
        for i in range( len( response ) ):
            if i < len( response ) - 1:
                response_diff_m.append( (response[i] + response[i+1]) / 2 )
            else:
                response_diff_m.append( 0 )

        l_response_int.append(response_diff_m)

    return l_response_int



# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
def write_response ( l_response, nombre ):

    encabezado = ''
    for elem in range( len( l_response ) ):
        if elem > 0 and elem < len( l_response ):
            encabezado += '|'
        encabezado += nombre + ' '
        encabezado +=  str( elem + 1 ) 

    try:
        np.savetxt( nombre + '.txt', np.array( l_response ).T, delimiter = "|", header = encabezado, comments = '', fmt='%.15f' )
        print( '\n\t [*] El archivo ' + nombre + '.txt' + ' se generó exitosamente. \n' )
    except:
        print ( '\n\t [!] Error al crear ' + nombre + '.txt \n' )



# > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
def main ():

    # acq_name = input( "\n\t Escribe el nombre del archivo: " )
    # acq_name = "R54 Ep + Lev SE 21 noviembre 2018 Registro 8 abril 2019 HD GD Señal 5 Curva IO Pre Lev.acq"
    
    formInicial = Tk()
    formInicial.title( "Analizador de señales" )
    formInicial.geometry("400x40")

    # Apertura del archivo acq
    acq_name = search_file()

    messagebox.showinfo( message = "La lectura del archivo podría tardar.", title = "Lectura ACQ" )   

    data = open_acq( acq_name )

    register = data.channels[0].data
    reference = data.channels[1].data
    
    response_size = 2000

    print(data.samples_per_second)
    
    response, l_response = get_response ( register, reference, response_size )
    l_response_fitted, response_fitted = filter( l_response, response )

    l_response_diff_r = diff_right( l_response_fitted )
    l_response_diff_m = diff_mean( l_response_fitted )
    l_response_int = int_trapezium( l_response_fitted )
    print ("\n\t  Se encontraron " + str(len(l_response)) + " respuestas.")

    write_response( l_response, "Respuesta" )
    write_response( l_response_fitted, "Respuesta_suave" )
    write_response( l_response_diff_r, "Diferencia_derecha" )
    write_response( l_response_diff_m, "Diferencia_central" )
    write_response( l_response_int, "Integral_Trapecio" )

    # plot_response( np.array(l_response).T, response, response_size, acq_name )
    plot_response_2( np.array(l_response).T, response, response_fitted, response_size, acq_name )

    formInicial.mainloop()
    
    return 0

main ()