# Программа расчёта 3-ёх мерного полёта ракеты.

# Процедура вывода данных о результате полёта в файл.
# Получает:		Данные о результате интегрирования полёта.
# Делает: 		Записывает в файл.
# Возвращает:		Ничего.
def answer_to_f(dget):
	pass	# Заглушка на время.

if __name__ == "__main__":
	# 1. Составные части.
	print("# 1. Составные части.")
	import windtower	# Содержит в себе функции полёта.
	import libf2		# Содержит в себе функции работы с файлами.
	#import 3d_main_h	# Содержит в себе функции подработки данных, чтобы удобно передать в целевую функцию.
	from sys import argv	# Нужно для получения данных при запуске.
	
	# 2. Получить данные argv.
	print("# 2. Распаковка данных argv.")
	script, first = argv	# Нас интересует раздел <first>. Через него будет передаваться флаг на исполнение программы.
	
	# 3. Загрузка config.
	print('# 3. Загрузка config.')
	dcfg = libf2.D_FROM_FILE('config')
	
	# 4. Загрузка данных.
	print("# 4. Загрузка данных.")
	ddata_rocket = libf2.D_FROM_FILE_FLOAT(dcfg['filename_rocket_input'])
	ddata_object = libf2.D_FROM_FILE_FLOAT(dcfg['filename_object_input'])
	
	# 5. Подготовка значений и запуск программы.
	print("# 5. Подготовка значений и запуск программы.")
	if(first == "to_point"):
		#dget = 3d_main_h.to_point(ddata_rocket, ddata_object)
		
		# Подготовка значений.
		dget = {}
		dget = ddata_rocket | ddata_object
		dret = windtower.rabor(dget['mu_0'], dget['mu_1'], dget['yakor_0'], dget['yakor_1'], dget)
		#print(dret)
		dget = dget | dret
		dget['if_povorot_start'] = dget['flag_startpovorot']
		dget['if_sbros_march'] = dget['flag_sbros_start']
		dget['dt'] = 0.01
		
		
		print("# ПУСК")
		# ПУСК.
		dfrom = windtower.d_to_point(dget, ddata_object['x_0'], ddata_object['y_0'], ddata_object['z_0'])
		
		print("# Анализ решения.")
		# Анализ решения.
		print(dfrom['answer'])
		print("x_r = %s,\ty_r = %s,\tz_r = %s" %(dfrom['l_x'][-1], dfrom['l_y'][-1], dfrom['l_z'][-1]))
		
		# 1. Расписать файл координат.
		f = open("DUMP_graphic_xyz.txt", "w")
		counter = 0
		for elem in dfrom['l_x']:
			f.write("%s;%s;%s\n" %(dfrom['l_z'][counter], elem, dfrom['l_y'][counter]))
			counter += 1
		f.close()
		
		# 2. Расписать файл масс.
		f = open("graphic_m.txt", 'w')
		counter = 0
		for elem in dfrom['l_t']:
			f.write("%s;%s\n" %(elem, dfrom['l_m'][counter]))
			counter += 1
		f.close()
		
		# 3. Расписать файлы скоростей.
		#print(len(dfrom['l_V_x']))
		#print(len(dfrom['l_V_y']))
		#print(len(dfrom['l_V_z']))
		#f = open("graphic_vx.txt", 'w')
		#counter = 0
		#for elem in dfrom['l_V_x']:
		#	f.write("%s;%s\n" %(dfrom['l_t'][counter], elem))
		#	counter += 1
		#f.close()
		#
		#counter = 0
		#f = open("graphic_vy.txt", 'w')
		#for elem in dfrom['l_V_y']:
		#	f.write("%s;%s\n" %(dfrom['l_t'][counter], elem))
		#	counter += 1
		#f.close()
		#
		#counter = 0
		#f = open("graphic_vz.txt", 'w')
		#for elem in dfrom['l_V_z']:
		#	f.write("%s;%s\n" %(dfrom['l_t'][counter], elem ))
		#	counter += 1
		#f.close()
		f = open("DUMP_graphic_vsum.txt", 'w')
		counter = 0
		for elem in dfrom['l_V_sum']:
			f.write("%s;%s\n" %(dfrom['l_t'][counter], elem))
			counter += 1
		f.close()
	elif(first == "dog"):
		letter = "dog"
		# Подготовка значений.
		dget = {}
		dget = ddata_rocket | ddata_object
		dret = windtower.rabor(dget['mu_0'], dget['mu_1'], dget['yakor_0'], dget['yakor_1'], dget)
		#print(dret)
		dget = dget | dret
		dget['if_povorot_start'] = dget['flag_startpovorot']
		dget['if_sbros_march'] = dget['flag_sbros_start']
		dget['dt'] = 0.01
		
		
		print("# ПУСК")
		# ПУСК.
		dfrom = windtower.d_to_ALL(letter, dget, ddata_object['x_0'], ddata_object['y_0'], ddata_object['z_0'])
		
		print("# Анализ решения.")
		# Анализ решения.
		print(dfrom['answer'])
		print("x_r = %s,\ty_r = %s,\tz_r = %s" %(dfrom['l_x'][-1], dfrom['l_y'][-1], dfrom['l_z'][-1]))
		
		# 1. Расписать файл координат.
		f = open("DUMP_graphic_xyz.txt", "w")
		counter = 0
		for elem in dfrom['l_x']:
			f.write("%s;%s;%s\n" %(dfrom['l_z'][counter], elem, dfrom['l_y'][counter]))
			counter += 1
		f.close()
		
		# 2. Расписать файл масс.
		f = open("graphic_m.txt", 'w')
		counter = 0
		for elem in dfrom['l_t']:
			f.write("%s;%s\n" %(elem, dfrom['l_m'][counter]))
			counter += 1
		f.close()
		
		# 3. Расписать файлы скоростей.
		#print(len(dfrom['l_V_x']))
		#print(len(dfrom['l_V_y']))
		#print(len(dfrom['l_V_z']))
		#f = open("graphic_vx.txt", 'w')
		#counter = 0
		#for elem in dfrom['l_V_x']:
		#	f.write("%s;%s\n" %(dfrom['l_t'][counter], elem))
		#	counter += 1
		#f.close()
		#
		#counter = 0
		#f = open("graphic_vy.txt", 'w')
		#for elem in dfrom['l_V_y']:
		#	f.write("%s;%s\n" %(dfrom['l_t'][counter], elem))
		#	counter += 1
		#f.close()
		#
		#counter = 0
		#f = open("graphic_vz.txt", 'w')
		#for elem in dfrom['l_V_z']:
		#	f.write("%s;%s\n" %(dfrom['l_t'][counter], elem ))
		#	counter += 1
		#f.close()
		f = open("DUMP_graphic_vsum.txt", 'w')
		counter = 0
		for elem in dfrom['l_V_sum']:
			f.write("%s;%s\n" %(dfrom['l_t'][counter], elem))
			counter += 1
		f.close()
	elif(first=='fields'):
		letter = "dog"
		# Подготовка значений.
		dget = {}
		dget = ddata_rocket | ddata_object
		dret = windtower.rabor(dget['mu_0'], dget['mu_1'], dget['yakor_0'], dget['yakor_1'], dget)
		#print(dret)
		dget = dget | dret
		dget['if_povorot_start'] = dget['flag_startpovorot']
		dget['if_sbros_march'] = dget['flag_sbros_start']
		dget['dt'] = 0.01
		
		# Создание полей.
		l_xcto = []	# Списки координат начальных цели.
		l_ycto = []
		l_zcto = []
		limx = dget['field_limx']
		limy = dget['field_limy']
		limz = dget['field_limz']
		iter_x = 0
		while(iter_x <= limx):
			iter_y = 0
			while(iter_y <= limy):
				iter_z = -limz
				while(iter_z <= limz):
					xc = iter_x
					yc = iter_y
					zc = iter_z
					print("# ПУСК по цели: xc = %s\tyc = %s\t zc = %s" %(xc, yc, zc))
					dfrom = windtower.d_to_ALL(letter, dget, xc, yc, zc)
					if(dfrom['answer']):	# Если вдруг решение прошло, то его следует сохранить.
						print(dfrom['answer'])
						l_xcto.append(iter_x)
						l_ycto.append(iter_y)
						l_zcto.append(iter_z)
					iter_z += 10000
				iter_y += 10000
			iter_x += 10000
		# Поворот задачи.
		#a = input("=============================================================")
		dget['v_0x'] = dget['v_0x'] * -1
		#print(dget['v_0x'])
		l_xcon = []
		l_ycon = []
		l_zcon = []
		iter_x = 0
		while(iter_x <= limx):
			iter_y = 0
			while(iter_y <= limy):
				iter_z = -limz
				while(iter_z <= limz):
					xc = iter_x
					yc = iter_y
					zc = iter_z
					print("# ПУСК по цели: xc = %s\tyc = %s\t zc = %s" %(xc, yc, zc))
					dfrom = windtower.d_to_ALL(letter, dget, xc, yc, zc)
					if(dfrom['answer']):	# Если вдруг решение прошло, то его следует сохранить.
						print(dfrom['answer'])
						l_xcon.append(iter_x*-1)
						l_ycon.append(iter_y*-1)
						l_zcon.append(iter_z*-1)
						#print(iter_x*-1)
					iter_z += 10000
				iter_y += 10000
			iter_x += 10000
		
		# Анализ решения.
		print("# Анализ решения.")
		
		# 1. Расписать файл координат.
		#f = open("DUMP_graphic_xyz.txt", "w")
		#counter = 0
		#for elem in l_xc:
		#	f.write("%s;%s;%s\n" %(dfrom['l_z'][counter], elem, dfrom['l_y'][counter]))
		#	counter += 1
		#f.close()
		
		# 2. Расписать файл масс.
		#f = open("graphic_m.txt", 'w')
		#counter = 0
		#for elem in dfrom['l_t']:
		#	f.write("%s;%s\n" %(elem, dfrom['l_m'][counter]))
		#	counter += 1
		#f.close()
		
		# 3. Расписать файлы проекций плоскостей.
		# 3.1 XY
		f = open("DUMP_graphic_field_xy.txt", 'w')
		counter = 0
		for elem in l_xcto:
			f.write("%s;%s\n" %(elem, l_ycto[counter]))
			counter += 1
		counter = 0
		for elem in l_xcon:
			f.write("%s;%s\n" %(elem, l_ycon[counter]))
			counter += 1
		f.close()
		
		# 3.2 XZ
		f = open("DUMP_graphic_field_xz.txt", 'w')
		counter = 0
		for elem in l_xcto:
			f.write("%s;%s\n" %(elem, l_zcto[counter]))
			counter += 1
		counter = 0
		for elem in l_xcon:
			f.write("%s;%s\n" %(elem, l_zcon[counter]))
		f.close()
		
		# 3.3 YZ
		f = open("DUMP_graphic_field_yz.txt", 'w')
		counter = 0
		for elem in l_ycto:
			f.write("%s;%s\n" %(elem, l_zcto[counter]))
			counter += 1
		counter = 0
		for elem in l_ycon:
			f.write("%s;%s\n" %(elem, l_zcon[counter]))
		f.close()
