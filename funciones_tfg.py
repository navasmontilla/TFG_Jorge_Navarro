def crear_mapa():
    import leafmap

    # Crear el mapa
    Mapa = leafmap.Map(center=(40.0, -3.5), zoom=6)

    # URL del WMS (sin SERVICE=WMS) - IGN elevaciones
    wms_url = "https://servicios.idee.es/wms-inspire/mdt?"

    # Listar capas disponibles (para comprobar si "EL.ElevationGridCoverage" es correcta)
    leafmap.get_wms_layers(wms_url)

    # Añadir capa WMS de elevaciones
    Mapa.add_wms_layer(
        url=wms_url,
        layers="EL.ElevationGridCoverage",
        name="Elevation Grid Coverage",
        opacity=1.0,
        shown=False
    )

    # URL del WMS - Forestal
    wms_url = "https://wms.mapama.gob.es/sig/Biodiversidad/MFE/wms.aspx?"

    # Listar capas disponibles (para comprobar si "LC.LandCoverSurfaces" es correcta)
    leafmap.get_wms_layers(wms_url)

    # Añadir capa WMS de cobertura del suelo
    Mapa.add_wms_layer(
        url=wms_url,
        layers="LC.LandCoverSurfaces",
        name="Land cover",
        opacity=1.0,
        shown=False
    )
    
    
    # URL del WMS - Ortofoto PNOA
    wms_url = "https://www.ign.es/wms-inspire/pnoa-ma"

    # Listar capas disponibles 
    leafmap.get_wms_layers(wms_url)

    # Añadir capa WMS de Ortoimagen
    Mapa.add_wms_layer(
        url=wms_url,
        layers="OI.OrthoimageCoverage",
        name="Ortoimage PNOA",
        opacity=1.0,
        shown=False
    )
    
    
    return Mapa

def extraer_caja_delimitadora(feature, base_res):
    import numpy as np
    coords = feature["geometry"]["coordinates"][0]
    lons = [c[0] for c in coords]
    lats = [c[1] for c in coords]

    min_lon, max_lon = min(lons), max(lons)
    min_lat, max_lat = min(lats), max(lats)
    aspectR = (max_lon - min_lon) / (max_lat - min_lat)
    
    center_lat = min_lat + (max_lat - min_lat) / 2.0
    f=  np.cos(np.radians(center_lat)) #factor corrector longitud: resoluciín mayor para longitud para tener malla cuadrada en distancias

    Nx = int((max_lon - min_lon) / (base_res/f))
    Ny = int((max_lat - min_lat) / base_res)

    max_lon = min_lon + Nx * (base_res/f)
    max_lat = min_lat + Ny * base_res 

    return min_lon, max_lon, min_lat, max_lat, Nx, Ny, aspectR
    
def generar_cuadricula(min_lon, max_lon, min_lat, max_lat, Nx, Ny):
    import numpy as np
    lons_array = np.linspace(min_lon, max_lon, Nx)
    lats_array = np.linspace(min_lat, max_lat, Ny)
    grid_points = [(lon, lat) for lat in lats_array for lon in lons_array]
    return grid_points
    
def obtener_elevaciones(grid_points, Nx, Ny):
    import json  # Para trabajar con datos JSON
    import requests  # Para enviar solicitudes HTTP a la API
    import numpy as np

    features = [{
        "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [lon, lat]}
    } for lon, lat in grid_points]

    geojson_payload = {
        "type": "FeatureCollection",
        "features": features
    }

    payload = {
        "inputs": {
            "crs": 4326,
            "formato": "geojson",
            "geom": json.dumps(geojson_payload),
            "outputFormat": "array",
            "withCoord": False
        }
    }

    url = "https://api-processes.idee.es/processes/getElevation/execution"
    response = requests.post(url, json=payload)

    if response.status_code == 200:
        data = response.json()
        return np.array(data["values"]).reshape(Ny, Nx)
    else:
        raise Exception("Error al obtener elevaciones")


def obtener_cobertura_vegetal(bbox, size, wms_url, wms_layer):
    from owslib.wms import WebMapService
    from rasterio.io import MemoryFile

   
    wms = WebMapService(wms_url)
    response = wms.getmap(
        layers=[wms_layer],
        srs="EPSG:4326",
        bbox=bbox,
        size=size,
        format="image/tiff",
        transparent=True
    )
    with MemoryFile(response.read()) as memfile:
        with memfile.open() as dataset:
            return dataset.read(1)


def procesar_dominio(features, base_res):
    import numpy as np 
    
    if not features:
        print("Aún no se han dibujado formas.")
        return

    feature = features[0]
    marker_feature = features[1] if len(features) > 1 else None

    if feature["geometry"]["type"] != "Polygon":
        print("Por favor, dibuja un rectángulo.")
        return

    min_lon, max_lon, min_lat, max_lat, Nx, Ny, aspectR = extraer_caja_delimitadora(feature, base_res)
    
    # Ahora vamos a generar una malla en km
    lon_grid = np.linspace(min_lon, max_lon, Nx)  # Longitudes desde xmin hasta xmax
    lat_grid = np.linspace(min_lat, max_lat, Ny)  # Latitudes desde ymin hasta ymax
    center_lon = min_lon + (max_lon - min_lon) / 2.0
    center_lat = min_lat + (max_lat - min_lat) / 2.0
    lon_distances = (lon_grid - min_lon) * 111 * np.cos(np.radians(center_lat))  # Convertir grados a km
    lat_distances = (lat_grid - min_lat) * 111  # Convertir grados a km
    
    
    print(f"Caja delimitadora: ({min_lon}, {min_lat}) hasta ({max_lon}, {max_lat})")
    print(f"Relación de aspecto: {aspectR}, Nx: {Nx}, Ny: {Ny}")        
    
    if marker_feature and marker_feature["geometry"]["type"] == "Point":
        marker_lon, marker_lat = marker_feature["geometry"]["coordinates"]
        print(f"Coordenadas del marcador: ({marker_lon}, {marker_lat})")

    return Nx, Ny, min_lon, max_lon, min_lat, max_lat, center_lat, marker_lon, marker_lat, lon_distances, lat_distances    
            

def procesar_topografia(fact,Nx,Ny,min_lon, max_lon, min_lat, max_lat, marker_lon, marker_lat):
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    import numpy as np

    #Factor de escala
    Nx_fact = np.round(Nx/fact).astype(int)
    Ny_fact = np.round(Ny/fact).astype(int)
    
    #Calculo de una matriz mas pequeña para evitar fallos en la representación 
    grid_points_fact = generar_cuadricula(min_lon, max_lon, min_lat, max_lat, Nx_fact, Ny_fact)
    elevaciones_fact = obtener_elevaciones(grid_points_fact, Nx_fact, Ny_fact)
   
    #Calculo de los puntos reales 
    grid_points = generar_cuadricula(min_lon, max_lon, min_lat, max_lat, Nx, Ny)
    
    #Aplano los valores factorizados para poder aplicar griddata de (Ny_fact, Nx_fact) → (Ny_fact * Nx_fact,) Lo aplano para que este en 1D y se pueda usar 
    elevaciones_fact_flat = elevaciones_fact.flatten()#Hay que aplanarlo 
    
    print(Nx_fact,Ny_fact)
    #print(elevaciones_fact_flat)
    #print(grid_points_fact)
    
    # Interpolación con griddata
    elevaciones_interp = griddata(points=grid_points_fact, values=elevaciones_fact_flat, xi=grid_points, method='cubic') #MIRAR CODIGO EXPLICADO EN EL CORRE (FECHA CORREO 20/05)
    
    # Convertimos el resultado en forma de matriz para visualización
    elevaciones = np.array(elevaciones_interp).reshape(Ny, Nx)

    return elevaciones
    
    
def graficar_elevacion(elevaciones,lon_distances,lat_distances):
    import numpy as np
    import matplotlib.pyplot as plt

    # Crear el gráfico usando Matplotlib
    plt.figure(figsize=(8, 6))  # Crear una nueva figura con un tamaño específico (8x6 pulgadas)

    # Representar los datos de elevación usando imshow, que mostrará las elevaciones como una imagen 2D
    # El argumento 'extent' define los límites del gráfico usando la caja delimitadora
    # 'cmap="terrain"' aplica un mapa de colores adecuado para elevaciones del terreno (e.g., elevaciones altas en marrón/amarillo)
    # 'origin="lower"' coloca el origen del gráfico en la esquina inferior izquierda (la latitud aumenta hacia arriba)
    plt.imshow(elevaciones, extent=[lon_distances[0], lon_distances[-1], lat_distances[0], lat_distances[-1]], cmap="terrain", origin="lower")

    # plt.contour(elevations, extent=[lon_distances[0], lon_distances[-1], lat_distances[0], lat_distances[-1]], levels=30, colors="k", origin="lower")

    # Añadir una barra de color al lado del gráfico para indicar los valores de elevación
    plt.colorbar(label="Elevación (m)")  # Etiquetar la barra de color como Elevación en metros

    # Etiquetar el eje x (Longitud en km) y el eje y (Latitud en km)
    plt.xlabel("Longitud (km)")
    plt.ylabel("Latitud (km)")

    # Añadir un título al gráfico
    plt.title("Mapa de Elevación")

    # Mostrar el gráfico en pantalla
    plt.show()
    
    return
    


def obtener_imagen_cobertura_vegetal(min_lon, max_lon, min_lat, max_lat, img_width, img_height, base_res):
    from owslib.wms import WebMapService
    from PIL import Image
    import io
    import numpy as np
    
    # Cobertura vegetal
    wms_url = "https://wms.mapama.gob.es/sig/Biodiversidad/MFE/wms.aspx?"  # Usar esta URL

    # Conectar al servicio WMS
    wms = WebMapService(wms_url, version='1.3.0')

    bbox = (min_lon, min_lat, max_lon, max_lat)
    wms_layer = "LC.LandCoverSurfaces"  # Usar esta capa

    fc = max(int(base_res / 0.0001), 1)  # Factor corrector para resolución

    # Obtener datos raster del WMS
    response = wms.getmap(
        layers=[wms_layer],
        srs="EPSG:4326",
        bbox=bbox,
        size=(fc * img_width, fc * img_height),
        format="image/png",
        transparent=True
    )

    # Convertir la imagen a un array NumPy
    image = Image.open(io.BytesIO(response.read()))
    image_np = np.array(image)

    return image_np, fc

def resize_and_classify(image_np, fc):
    from owslib.wms import WebMapService
    from PIL import Image
    import io
    import numpy as np
    
    # Reescalar la matriz genérica (x, y, 4) a (x, y, 1)
    if image_np.shape[2] != 4:
        raise ValueError("La imagen debe tener 4 canales (RGBA).")

    # Reescalar la imagen manteniendo la relación de aspecto
    resized_image = image_np[::fc, ::fc]

    # Clasificación de píxeles según los valores de color
    color_map = {
        (85, 185, 51, 255): 1,   # Arbolado
        (187, 214, 132, 255): 2,  # Arbolado ralo
        (143, 143, 51, 255): 3,   # Arbolado disperso
        (238, 255, 203, 255): 4,  # Desarbolado
        (255, 145, 0, 255):5,     # Humedades
        (255, 236, 181, 255): 6,  # Cultivos
        (194, 194, 194, 255): 7,  # Artificial
        (143, 243, 248, 255): 8   # Agua
    }

    rows, cols, channels = resized_image.shape
    classified_matrix = np.zeros((rows, cols), dtype=int)

    for i in range(rows):
        for j in range(cols):
            pixel = tuple(resized_image[i, j])  # Tomar el valor de cada píxel
            classified_matrix[i, j] = color_map.get(pixel, 0)  # Asignar valor o 0 si no está definido

    return resized_image, classified_matrix



def mostrar_cobertura_con_leyenda(resized_image):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    
    # Mostrar correspondencia entre valores y tipos de cobertura vegetal
    legend_labels = {
        "Arbolado": (85, 185, 51),
        "Arbolado ralo": (187, 214, 132),
        "Arbolado disperso": (143, 143, 51),
        "Desarbolado": (238, 255, 203),
        "Humedades": (255, 145, 0),
        "Cultivos": (255, 236, 181),
        "Artificial": (194, 194, 194),
        "Agua": (143, 243, 248)
    }

    # Crear la leyenda
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.imshow(resized_image)
    ax.axis("off")  # Ocultar ejes
    ax.set_title("Cobertura Vegetal con Leyenda")

    # Añadir la leyenda
    patches = [mpatches.Patch(color=np.array(color) / 255, label=label) for label, color in legend_labels.items()]
    plt.legend(handles=patches, loc="lower right", title="Leyenda")
    plt.show()


def procesar_cobertura_vegetal(min_lon, max_lon, min_lat, max_lat, img_width, img_height, base_res):
    import numpy as np
    
    # Paso 1: Obtener imagen y factor corrector
    image_np, fc = obtener_imagen_cobertura_vegetal(min_lon, max_lon, min_lat, max_lat, img_width, img_height, base_res)

    # Paso 2: Clasificar
    resized_image, result = resize_and_classify(image_np, fc)

    # Paso 3: Mostrar
    mostrar_cobertura_con_leyenda(resized_image)

    

    result = np.flipud(result)

    return result

    
    
def crear_matriz_xx(result, valor_xx):
    import numpy as np

    matriz = np.array(result, dtype=np.float32)

    matriz_xx = np.zeros_like(matriz)

    # Asignar valores usando el vector recibido
    matriz_xx[matriz == 1] = valor_xx[0]  # Arbolado
    matriz_xx[matriz == 2] = valor_xx[1]  # Arbolado disperso
    matriz_xx[matriz == 3] = valor_xx[2]  # Arbolado ralo
    matriz_xx[matriz == 4] = valor_xx[3]  # Desarbolado
    matriz_xx[matriz == 5] = valor_xx[4]  # Humedades
    matriz_xx[matriz == 6] = valor_xx[5]  # Cultivos
    matriz_xx[matriz == 7] = valor_xx[6]  # Artificial
    matriz_xx[matriz == 8] = valor_xx[7]  # Agua

    return matriz_xx

    

def generar_input_simulacion(lon_distances, lat_distances, Nx, Ny, elevaciones, marker_lon, marker_lat, matriz_sc, min_lon, max_lon, min_lat, max_lat, center_lat,radio,u_x,u_y, K, h, Tinf, RH,L, matriz_rfmax,matriz_rhof,matriz_cf,matriz_h,matriz_mg,matriz_alpha,matriz_tpc,matriz_a,FinalTime,DumpTime):
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    
    def configuracion_rutas():
        #ADAPTACION DE FIRESET UP A COLAB
        global folder_case, folder_out
        global fname_config, fname_configF, fname_initialF, fname_landscapeF

        folder_case="case"
        folder_out="output-files"

        fname_config="case/configure.input"
        fname_configF="case/configureFire.input"
        fname_initialF="case/initialFire.out"
        fname_landscapeF="case/landscapeFire.out"

    def configuracion_global():
        #Simulation setup - Global settings
        global SizeX, SizeY, SizeZ #FinalTime
        global CFL, Order
        global xcells, ycells
        global Face_1, Face_2, Face_3, Face_4
        global u_x, u_y, K, h, Tinf, RH, L

        #FinalTime = 2600.0  #this is the final time before 800
        SizeX = lon_distances[-1]*1000      #this is the domain size in X METROS
        SizeY = lat_distances[-1]*1000     #this is the domain size in Y METROS
        SizeZ = 100.0   ##??

        #Fire model setup
        CFL = 0.4          #Courant-Friedrichs-Lewy stability criterion. Must be <0.5 NO TOCAR
        Order = 7          #Order of accuracy for advection. Should be 5 or 7.

        #Mesh setup
        xcells = Nx       #number of cells in X
        ycells = Ny       #number of cells in Y

        #Boundary conditions (do not modify them!)
        Face_1 = 3 #-y
        Face_2 = 3 #+x
        Face_3 = 3 #+y
        Face_4 = 3 #-x



    def inicializar_mallas_y_arrays():
        global dx, dy, T, Y, rhof, Cf, H, Z, alpha, Tpc, SC, Rfmax, Mg, A
        global x, y, xc, yc

        dx=SizeX/xcells
        dy=SizeY/ycells

        T=np.zeros((xcells,ycells))
        Y=np.zeros((xcells,ycells))
        rhof=np.zeros((xcells,ycells))
        Cf=np.zeros((xcells,ycells))
        H=np.zeros((xcells,ycells))
        Z=np.zeros((xcells,ycells))
        alpha=np.zeros((xcells,ycells))
        Tpc=np.zeros((xcells,ycells))
        SC=np.zeros((xcells,ycells))
        Rfmax=np.zeros((xcells,ycells))
        Mg=np.zeros((xcells,ycells))
        A=np.zeros((xcells,ycells))

        x=np.arange(0+dx/2.0, SizeX, dx)
        y=np.arange(0+dy/2.0, SizeY, dy)

        xc, yc= np.meshgrid(x,y,indexing='ij')

        print(Z.shape,elevaciones.shape)
        print(SC.shape, matriz_sc.shape)
        

    def establecer_condiciones_iniciales():
        marker_lon_km = (marker_lon - min_lon) * 111 * np.cos(np.radians(center_lat))  # Convertir grados a km
        marker_lat_km = (marker_lat - min_lat) * 111  # Convertir grados a km
        
       
        

        xcenter = marker_lon_km*1000 #convertir a metros
        ycenter= marker_lat_km*1000 #convertir a metros

        for l in range(0,xcells):
            for m in range(0,ycells):
                r=np.sqrt((xc[l,m]-xcenter)**2+(yc[l,m]-ycenter)**2)
                if r<radio:
                    T[l,m]=670.0
                else:
                    T[l,m]=300.0 #POner 300, es una prueba

                Y[l,m]=1.0
                
                #DEFINIR UNA MATRIZ PARA CADA ELEMENTO EN COVERTURA VEGETAL 
                Rfmax[l,m]=matriz_rfmax[m,l] #maximum volumetric porosity RCDC(relacionados cantidad de combustible), MODIFICABLE
                SC[l,m]= matriz_sc[m,l] #surface coverage # 
                Z[l,m]= elevaciones[m,l]   #TOPOGRAFIA 11/04
                rhof[l,m]=matriz_rhof[m,l]  #solid fuel density
                Cf[l,m]=matriz_cf[m,l]   #solid fuel specific heat NO TOCAR
                H[l,m]=matriz_h[m,l]    #solid fuel heat power PODER CALORIFICO DEL COMBUSTIBLE, ENERGIA QUE SE LIBEROA POR CADA KG DE COMBUSTIBLE, MODIFICABLE
                Mg[l,m]=matriz_mg[m,l]    #green wood water content ENTRE 0 Y 2, MODIFICABLE
                alpha[l,m]=matriz_alpha[m,l]  #ratio between green and season wood ALPHA 1, TODO VIVO.. ALPHA 0 , TODO MUERTO 
                Tpc[l,m]=matriz_tpc[m,l]   #pirolisis temperature TEMPERATURA A LA Q SE ACTIVA EL FUEGO, MODIFICABLE
                A[l,m]=matriz_a[m,l]   #reaction rate MUY SENSIBLE, MODIFICABLE


    def escribir_configuracion_y_salidas():
        f = open(fname_config, "w")
        f.write("/////SIMULATION_SETUP////// \n")
        f.write("FinalTime    "+str(FinalTime)+"\n")
        f.write("SizeX    "+str(SizeX)+"\n")
        f.write("SizeY    "+str(SizeY)+"\n")
        f.write("SizeZ    "+str(SizeZ)+"\n")
        f.close()

        f = open(fname_configF, "w")
        f.write("/////SIMULATION_SETUP////// \n")
        f.write("DumpTime    "+str(DumpTime)+"\n")
        f.write("CFL    "+str(CFL)+"\n")
        f.write("Order    "+str(Order)+"\n")
        f.write("////////MESH_SETUP/////////\n")
        f.write("xcells    "+str(xcells)+"\n")
        f.write("ycells    "+str(ycells)+"\n")
        f.write("///////BOUNDARY_COND///////\n")
        f.write("Face_1    "+str(Face_1)+"\n")
        f.write("Face_2    "+str(Face_2)+"\n")
        f.write("Face_3    "+str(Face_3)+"\n")
        f.write("Face_4    "+str(Face_4)+"\n")
        f.write("/////////PARAMETERS////////\n")
        f.write("u_x(m/s)    "+str(u_x)+"\n")
        f.write("u_y(m/s)    "+str(u_y)+"\n")
        f.write("rho(kg/m^3)    40.0 \n")
        f.write("C(kJ/(kg·K))   1.0 \n")
        f.write("H(kJ/kgFuel)   4000.0\n")
        f.write("K(kW/(m·K))    "+str(K)+"\n")
        f.write("h(kW/(m^3·K))    "+str(h)+"\n")
        f.write("Tinf(K)    "+str(Tinf)+"\n")
        f.write("Tpc(K)    400\n")
        f.write("RH    "+str(RH)+"\n")
        f.write("L(m)    "+str(L)+"\n")
        f.close()

        f = open(fname_initialF, "w")
        f.write("VARIABLES = X, Y, T, Y \n")
        f.write("CELLS = "+str(xcells)+", "+str(ycells)+"\n")
        for l in range(0,xcells):
            for m in range(0,ycells):
                f.write(str(xc[l,m])+" "+str(yc[l,m])+" "+str(T[l,m])+" "+str(Y[l,m])+"\n")
        f.close()

        f = open(fname_landscapeF, "w")
        f.write("VARIABLES = X, Y, Z, Rfmax, SC, rhof, Cf, Mg, alpha, Tp, H, A \n")
        f.write("CELLS = "+str(xcells)+", "+str(ycells)+"\n")
        for l in range(0,xcells):
            for m in range(0,ycells):
                f.write(str(xc[l,m])+" "+str(yc[l,m])+" "+str(Z[l,m])+" "+str(Rfmax[l,m])+" "+str(SC[l,m])+" "+str(rhof[l,m])+" "+str(Cf[l,m])+" "+str(Mg[l,m])+" "+str(alpha[l,m])+" "+str(Tpc[l,m])+" "+str(H[l,m])+" "+str(A[l,m])+"\n")
        f.close()

    def graficar_variables():
        print(Z.shape, SC.shape)
        
        #SC_rotada =  np.fliplr(SC) #Efecto espejo por el eje X, hay que hacerlo porque como esta programado guarda los ejes mal, al girar 90º a izquierda se orientan bien 
        
        plt.figure(figsize=(20, 8))
        plt.suptitle("Fire Model Variables", fontsize=16)

        #ORIGINAL variables = [T, Y, rhof, Cf, H, Z, alpha, Tpc, SC, Rfmax]
        variables = [T, Y, rhof, Cf, H, Z, alpha, Tpc, SC, Rfmax]
        titles = ["Temperature (K)", "Biomass Fraction", "Fuel Density (kg/m³)", "Specific Heat",
                  "Heat Power (kJ/kg)", "Topography (m)", "Green to dead wood ratio", "Pyrolysis Temperature (K)",
                  "Surface Coverage", "Rf0"]

        for i in range(10):
            
            ax = plt.subplot(2, 5, i + 1)
            contour = ax.contourf(xc, yc, variables[i], cmap="viridis")
            ax.set_aspect('equal')
            ax.set_title(titles[i])
            ax.set_xlabel("X (m)")
            ax.set_ylabel("Y (m)")
            plt.colorbar(contour, shrink=0.8)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig("fire_model_variables.png", dpi=300)

        print("Files written")

    # Llamada a todas las funciones internas
    configuracion_rutas()
    configuracion_global()
    inicializar_mallas_y_arrays()
    establecer_condiciones_iniciales()
    escribir_configuracion_y_salidas()
    graficar_variables()
