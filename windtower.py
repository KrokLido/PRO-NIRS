#	Файл с процедурами обсчёта пространственного полёта.
#	
#	Содержит:
#
#	<bellkeeper()>		# Процедура вычисления характеристик на каждом конкретном шаге.
#	<rabor()>		# Процедура вычисления геометрии аппарата.
#	<3d_to_point()>		# Процедура интегрирования при стрельбе ракетой в точку с условием того, что та должна повернуть.
#					точка полётная задаётся снаружи и передаётся внутрь функции.
#	Зависимости от своих модулей:
#				#	cxyz.py		- Расчёт коэфф-ов Сх, Су, Сz
#				#	graphic.py	- Строительство графиков.





#	Процедура осуществляет первоначальный расчёт массово-габаритных харк-к ЛА.
#	Получает:		- mu_0	-	относительная масса топлива стартовика
#				- mu_1	-	-//- маршевой ступени
#				- yakor_0	-	тяговооружённость стартовика
#				- yakor_1	-	тяговооружённость маршевой ступени
def rabor(mu_0, mu_1, yakor_0, yakor_1, ddata):
	d = ddata['d']
	F_midelya = 3.14 * d**2 /4
	
	betta = ddata['betta_0']
	g = 9.86
	
	m_sum = (ddata['m_payload']) / ((1 - betta*mu_0) * (1 - betta * mu_1) )
	
	m_0 = m_sum * mu_0 * betta
	m_0fuel = m_sum * mu_0
	
	m_1 = (m_sum - m_0 - ddata['m_payload'])
	m_1fuel = (m_sum - m_0) * mu_1

	P_10 = yakor_0 * m_sum * g
	P_11 = yakor_1 * (m_1 + ddata['m_payload']) * g

	G_10 = P_10 / ddata['I_10']
	G_11 = P_11 / ddata['I_11']

	ddata['mu_0'] = mu_0
	ddata['mu_1'] = mu_1
	ddata['F_midelya'] = F_midelya
	ddata['m_sum'] = m_sum
	ddata['m_0'] = m_0
	ddata['m_0fuel'] = m_0fuel
	ddata['m_1'] = m_1
	ddata['m_1fuel'] = m_1fuel
	ddata['P_10'] = P_10
	ddata['P_11'] = P_11
	ddata['G_10'] = G_10
	ddata['G_11'] = G_11
	ddata['t_0'] = m_0fuel/G_10
	ddata['t_1'] = m_1fuel/G_11
	return(ddata)





#	Процедура содержит в себе все формулы для вычисления всех координат по временному промежутку. 
#	Процедура вычисляет полёт с креном и скольжением.
#	Получает:		1. Данные с прошлого промежутка в виде словарика:
#					1.1 V_xy	-	Скорость поступательного движения x по xy
#					1.2 V_yy	-	Скорость поступательного движения y по xy
#					1.3 V_zx	-	Скорость поступательного движения х по zx
#					1.4 V_zz	-	Скорость поступательного движения z по zx
#					1.5 omega_x	-	Угловая скорость вращения вдоль оси x
#					1.6 omega_y	-	Угловая скорость вращения вдоль оси у
#					1.7 omega_z	-	Угловая скорость вращения вдоль оси z
#					1.8 gamma	-	Угол крена, т.е. угол вращения вокруг OX.
#					1.9 teta_sumxy	-	Суммарный угол скольжения ЛА в плоскости xy
#					1.10 teta_sumxz	-	Суммарный угол скольжения ЛА в плоскости xz
#					1.11 teta_xy	-	Угол атаки в плоскости OXY
#					1.12 teta_xz	-	Угол аатки в плоскости OXZ
#					1.13 atm	-	атмосфера
#					1.14 x		-	Координата x
#					1.15 y		-	Координата y
#					1.16 z		-	Координата z
#					1.17 rho	-	плотность атмосферы.
#					1.18 dt		-	временной промежуток.
#					1.19 V_ploskxy	-	Плоскостная скорость в OXY
#					1.20 V_ploskxz	-	Плоскостная скорость в OXZ
#					1.21 X_xy	-	АД-сила X_xy
#					1.22 Y_xy	-	АД-сила Y_xy
#					1.23 X_xz	-	АД-сила X_xz
#					1.24 Z_xz	-	АД-сила Z_xz
#					1.25 P_xy	-	Плоскостная тяга в OXY
#					1.26 P_xz	-	Плоскостная тяга в OXZ
#					1.27 m		-	Масса ЛА.
#					1.28 P		-	Переменная тяги. В зависимости от настроечного флага будет применяться P_10 else P_11
#						P_10 - тяга на старте.
#						P_11 - тяга на марше.
#					1.29 m_fuel	-	Масса топлива.
#					1.30 G		-	Переменная массового расхода. Зависимость от настроечного флага:
#						G_10 - расход на стартовом.
#						G_11 - расход на маршевом.
#				2. Таблицу данных по АД-коэффициентам. Более подробно смотреть в <cxyz.py>
#				3. Плотность воздуха на заданной высоте.
#				4. Настроечные флаги:
#					2.1 Пока без флагов.
#				5. Настроечный флаг применения массового расхода.
#					5.1 <False> - G_0 else G_1
#				6. Настроечный флаг выбора тяги.
#					6.1 <False> - P_0 else P_1
#	Делает:			Проводит вычисление.
#	Возвращает:		Результаты вычисления в виде словарика.
def bellkeeper(ddata, ldc, rho,  dt=0.01, flag_G=False, flag_P=False):
	from trigonometria import sin, cos, tg, atg
	import atmosphere_GOST4401_81 as atm
	import cxyz
	
	
	# Процедура осуществляет получение значения коэффициента.
	# Получает:			1. имя коэффициента
	#				2. таблицу всех коэффициентов.
	#				3. координату поиска. В НАШЕМ СЛУЧАЕ КООРДИНАТОЙ ЯВЛЯЕЕТСЯ УГОЛ!!!!
	# Делает:			Осуществляет выборку данных.
	# Возвращает:			Значение коэффициента.
	def xyzk(name, ldc, xy):
		if(xy == 0):
			return(0)
		elif(xy < 0):
			xy = abs(xy)
		if(name == "x"):
			lx=[]
			ly=[]
			lx = ldc[0]['lx'].copy()
			ly = ldc[0]['ly'].copy()
			table = {}
			table['x'] = lx
			table['y'] = ly
			cx = cxyz.get_value(table, xy)
			return(cx)
		elif(name == "y"):
			table = {}
			table['x'] =  ldc[1]['lx'].copy()
			table['y'] =  ldc[1]['ly'].copy()
			cy = cxyz.get_value(table, xy)
			return(cy)
		elif(name == "z"):
			lx=[]
			ly=[]
			table = {}
			table['x'] = ldc[2]['lx'].copy()
			table['y'] = ldc[2]['ly'].copy()
			cz = cxyz.get_value(table, xy)
			return(cz)
	
	# Получает:			1. Имя значения АД-силы.
	#				2. Таблицу АД-значений коэффициентов в заивисмости от некой координаты.
	#				3. Саму координату. В рассматриваемом случае координатой является угол скольжения.
	#				4. Плотность.
	#				5. Скорость подачи по данному направлению.
	# Делает:			Осуществляет расчёт значения АД-силы на данном направлении.
	# Возвращает:			Значение составляющей.
	def xyzforce(name, ldc, xy, rho, V):
		def force(cx, rho, v):
			return(cx * rho * v**2 / 2)
		
		if((name == 'x')or(name == "X")):
			cx = xyzk('x', ldc, xy)
			return(force(cx, rho, V))
		if((name == 'y')or(name == "Y")):
			cy = xyzk('y', ldc, xy)
			return(force(cy, rho, V))
		if((name == 'z')or(name == "Z")):
			cz = xyzk('y', ldc, xy)
			return(force(cz, rho, V))
	
	
	# ЗАТЫЧКА
	# Получает:			1. X_xy	-	лобовое сопротивление в OXY
	#				2. Y_xy	-	подъёмную силу в OXY
	#				3. X_xz	-	лобовое сопротивление в OXZ
	#				4. Z_xz	-	сносящую силу в OXZ
	#				5. teta_xy	-	угол скольжения в OXY
	#				6. teta_xz	-	угол скольжения в OXZ
	# Делает:			Вычисляет полное АД-сопротивление ракете.
	# Возвращает:			АД-сопротивление(магия, да?)
	def airresist(X_xy, Y_xy, X_xz, Z_xz, teta_xy, teta_xz):
		full_airresist = None
		return(full_airresist)
	
	# Выбор массового расхода. Я лентяй, поэтому тут.
	if(flag_G):
		ddata['G']  = ddata['G_11']
	else:
		ddata['G'] = ddata['G_10']
	
	if(flag_P):
		ddata['P'] = ddata['P_11']
		#print('P = ', ddata['P'])
	else:
		ddata['P'] = ddata['P_10']
		#print('P = ', ddata['P'])
	
	# Получает:			1. V_x?
	#				2. V_y? or V_z?
	# Делает:			Вычисляет суммарную плоскостную скорость.
	# Возвращает:			Суммарную плоскостную скорость.
	def V_plosk(V_x, V_y):
		return (V_x**2 + V_y**2)**0.5
	
	# Получает:			1. V_ploskxy
	#				2. V_ploskxz
	#				3. teta_sumxy
	#				4. teta_xy
	#				5. teta_sumxz
	#				6. teta_xz
	# Делает:		Вычисляет полную скорость.
	def V_sum(V_ploskxy, V_ploskxz, teta_sumxy, teta_xy, teta_sumxz, teta_xz):
		ret = ( (V_ploskxy*cos(teta_sumxy - teta_xy) + V_ploskxz*cos(teta_sumxz - teta_xz))**2 + (V_ploskxy*sin(teta_sumxy - teta_xy))**2 + (V_ploskxz*sin(teta_sumxz - teta_xz))**2 )**0.5
		return(ret)
	
	
	
	# Как выяснилось, тяга у ЛА одна на все три оси, поэтому и тяга должна быть одна.
	# Вопрос этот решается путём использования одной тяги и двух тригонометр.функций.
	# Поэтому P_xy, P_xz отныне используются в терминах P.
	
	# Шаг 0. Вычисление АД-сопротивления.
	ddata['X'] = 1.7*xyzforce("x", ldc, (ddata['nu_p']-ddata['nu_v']), ddata['rho'], ddata['V_sum'])
	ddata['Y'] = xyzforce("y", ldc, (ddata['nu_p']-ddata['nu_v']), ddata['rho'], ddata['V_sum'])
	
	ddata['V_x'] += (dt / ddata['m'])*(ddata['P']*cos(ddata['nu_p']) - ddata['X']*cos(ddata['nu_v']) - ddata['Y']*sin(ddata['nu_v']))
	ddata['V_y'] += (dt / ddata['m'])*(ddata['P']*sin(ddata['nu_p']) - ddata['X']*sin(ddata['nu_v']) + ddata['Y']*cos(ddata['nu_v']) - ddata['m']*ddata['glob_g'])
	#ddata['nu_v'] = atg(ddata['V_y'] / ddata['V_x'])
	ddata['V_sum'] = (ddata['V_x']**2 + ddata['V_y']**2)**0.5
	
	# Пересложение скоростей.
	
	
	# Координаты.
	ddata['x'] += ddata['V_x']*dt*cos(ddata['psi'])
	ddata['y'] += ddata['V_y']*dt
	ddata['z'] += ddata['V_x']*dt*sin(ddata['psi'])

	# Масса.
	ddata['m'] -= ddata['G'] * dt
	ddata['m_fuel'] -= ddata['G'] * dt

	return(ddata)




#	==	НА ТЕКУЩИЙ МОМЕНТ ЯВЛЯЕТ СОБОЙ ТЕСТОВЫЙ ОТСТРЕЛ РАКЕТЫ С ЦЕЛЬЮ ПРОВЕРКИ РАБОТЫ СИСТЕМЫ УРАВНЕНИЙ	==
#	Общие расчёты по пространственной задаче при полёте в точку.
#	Получает:			1) ddata - словарь данных,
#						2) x, y, z - координаты цели.
# 						3) flag_graph - Флаг формы возрата решения.
# 						4) flag_ret_all - Флаг возврата данных решения.
#	Делает:				Интегрирует.
#	Возвращает:			1)	Если флаг <flag_graph> опущен:
# 								Печать в консоль, сбор данных в контейнер и его возврат.
# 							Если поднят:
# 								Строит графики.
# 						2)	Если флаг <flag_ret_all> опущен:
# 								Возврат <False> если решение неудовлетворяет.
# 								Возврат <True> если решение удовлетворительно.
# 							Если флаг <flag_ret_all> поднят:
# 								Возврат словаря из вычисленных списков. 
#	Описание:	Процедура должна получить уже сформированный словарь характеристик ddata. Процедура возвращает
#	обратно не словарь <ddata>, а словарь <dret>, который содержит:
# 										- l_x = [] - список координат горизонтали.
# 										- l_y = [] - список координат вертикали.
# 										- l_t = [] - список координат времени.
# 										- l_V = [] - список скоростей.
# 										- l_m = [] - список масс.
# 										- l_resist = [] - список сопротивлений.
# 										- answer - строка состояния решения. Доступны следующие состояния:
# 												"Ракета достигла цели на стартовом участке."
# 												"Ракета достигла цели на маршевом участке."
# 												"Ракета достигла цели на пассивном участке."
# 												"Ракета не долетела."
def d_to_point(ddata, xc, yc, zc, flag_graph=False, flag_ret_all=False):
	from trigonometria import sin, cos
	import atmosphere_GOST4401_81 as atm
	import cxyz 
	def inf_dreturn(l_x, l_y, l_z, l_t, l_V_xy, l_V_yy, l_V_zx, l_V_zz, l_V_sum, l_m, answer, flag):
		dreturn = {}
		dreturn['l_x'] = l_x.copy()
		dreturn['l_y'] = l_y.copy()
		dreturn['l_z'] = l_z.copy()
		dreturn['l_t'] = l_t.copy()
		dreturn['l_V_xy'] = l_V_xy.copy()
		dreturn['l_V_yy'] = l_V_yy.copy()
		dreturn['l_V_zx'] = l_V_zx.copy()
		dreturn['l_V_zz'] = l_V_zz.copy()
		dreturn['l_V_sum'] = l_V_sum.copy()
		dreturn['l_m'] = l_m.copy()
		dreturn['answer'] = answer
		dreturn['flag'] = flag
		return(dreturn)
	
	# Формовка строки характеризующей решение:
	def inf_answer(position, retflag=False):
		if(position=='start'):
			if(retflag):
				return(True)
			return("Ракета достигла своей цели на стартовом участке.")
		if(position=='march'):
			if(retflag):
				return(True)
			return('Ракета достигла своей цели на маршевом участке.')
		if(position=='fatigue'):
			if(retflag):
				return(True)
			return('Ракета достигла своей цели на пассивном участке.')
		if(position=='fall'):
			if(retflag):
				return(False)
			return('Ракета не долетела.')
		if(position=="M<V"):
			if(retflag):
				return(False)
			return('Ракета не набрала нужной скорости.')
		if(position=="teta"):
			if(retflag):
				return(False)
			return("Ракета начала вести угол атаки более 10 градусов, что не есть хорошо.")
		if(position=='atmosphere'):
			if(retflag):
				return(False)
			return("Ракета вылетела за пределы атмосферы.")
		else:
			if(retflag):
				return(False)
			return('Здесь был бред')
		
	# Функция расчёта дальности до цели.
	def L_to_point(x, y, z, xp, yp, zp):
		xt = abs(xp - x)
		yt = abs(yp - y)
		zt = abs(zp - z)
		print('\txt = ',xt)
		print('\tyt = ',yt)
		print('\tzt = ',zt)
		L = (xt**2 + yt**2 + zt**2)**0.5
		return(L)
	
    	# Перед решением получаем начальные данные.
	print("WNDTWR= получение начальных данных.")
	#dra = rho()
	dra = atm.get_table(fcsv=True)	#	Таблица плотностей атмосферы в зависимости от высоты.
	altmax = atm.getmax()	#	Получили таблицу зависимости атмосферы.
	dcx = {}	# 	Получили таблицу зависимости коэффициентов Cx, Cy, Cz.
	dcy = {}
	dcz = {}
	dcx['lx'], dcx['ly'] = cxyz.get_table("cx")
	dcy['lx'], dcy['ly'] = cxyz.get_table("cy")
	dcz['lx'], dcz['ly'] = cxyz.get_table("cz")
	ldc = []
	ldc.append(dcx)
	ldc.append(dcy)
	ldc.append(dcz)
	
	
    	
	# Перед решением снимаем копию данных для удобства обращения.
	#print(ddata)
	mu_0 = ddata['mu_0']
	mu_1 = ddata['mu_1']
	d = ddata['d']
	F_midelya = ddata['F_midelya']
	m_sum = ddata['m_sum']
	m_0 = ddata['m_0']
	m_0fuel = ddata['m_0fuel']
	m_1 = ddata['m_1']
	m_1fuel = ddata['m_1fuel']
	P_10 = ddata['P_10']
	P_11 = ddata['P_11']
	G_10 = ddata['G_10']
	G_11 = ddata['G_11']
	limpsi = ddata['limpsi']
	limnu = ddata['limnu']
	limgamma = ddata['limgamma']
	
	# Принимаем начальные данные.
	dt = ddata['dt']		# Шаг интегрирования по времени.
	t = 0			# Общее время.
	epsilon_m = 5	# Ошибка определения масс топлива.
	x = 0			# Начальная горизонтальная координата <x>.
	y = 0			# Начальная вертикальная координата.
	z = 0			# Начальная координата горизонта <z>.
	V = 0			# Начальная скорость.
	epsilon = 30	# Вероятное отклонение от цели, засчитываемое за попадание.
    	
	l_m = []
	l_vxy = []
	l_vyy = []
	l_vzx = []
	l_vzz = []
	l_vsum = []
	l_t = []
	l_x = []
	l_y = []
	l_z = []
	l_m.append(m_sum)
	l_vxy.append(0)
	l_vyy.append(0)
	l_vzx.append(0)
	l_vzz.append(0)
	l_vsum.append(0)
	l_t.append(0)
	l_x.append(0)
	l_y.append(0)
	l_z.append(0)
	
	# Теперь нужно сформировать словарь, который будет гоняться в функции <bellkeeper>.
	dbell = {}
	
	# Скорости.
	#dbell['V_xy'] = 0	# Скорости поступательного движения по осям.
	#dbell['V_yy'] = 0
	#dbell['V_zx'] = 0
	#dbell['V_zz'] = 0
	dbell['V_x'] = 0
	dbell['V_y'] = 0
	dbell['V_sum'] = 0
	
	# Массы.
	dbell['m'] = m_sum	# Масса
	dbell['m_fuel'] = m_0fuel + m_1fuel
	
	# Тяги.
	dbell['P'] = P_10	# Значения тяги.
	dbell['P_10'] = P_10
	dbell['P_11'] = P_11
	dbell['G_10'] = G_10
	dbell['G_11'] = G_11
	
	# Углы полёта.
	dbell['teta_sumxy'] = ddata['teta_0']	# Суммарный угол старта ракеты в вертикальной плоскости.
	dbell['teta_sumxz'] = ddata['psi']		# Суммарный угол старта ракеты в горизонтальной плоскости.
	dbell['teta_xy'] = 3	# Значение угла скольжения ЛА. Вертикальная плоскость.
	dbell['teta_xz'] = 0	# Значения угла скольжения ЛА в горизонтальной плоскости.
	dbell['omega_x'] = 0	# Угловая скорость.
	dbell['omega_y'] = 0	
	dbell['omega_z'] = 0
	dbell['gamma'] = ddata['gamma']		# Угол крена стартовый, обычно, равен нулю.
	
	dbell['nu_p'] = ddata['teta_0']
	dbell['nu_v'] = ddata['teta_0']
	dbell['psi'] = ddata['psi']
	
	# Значения АД-сил.
	dbell['X_xy'] = 0		# Значение АД-сил.
	dbell['Y_xy'] = 0
	dbell['X_xz'] = 0
	dbell['Z_xz'] = 0
	dbell['X'] = 0
	dbell['Y'] = 0
	
	# Значения моментов инерции и прочего.
	dbell['I_x'] = ddata['I_x']	# Моменты инерции вдоль осей движения.
	dbell['I_y'] = ddata['I_y']
	dbell['I_z'] = ddata['I_z']
	dbell['M_x'] = ddata['M_x']	# Балансировочные зависимости ЛА.
	dbell['M_y'] = ddata['M_y']
	dbell['M_z'] = ddata['M_z']
	
	dbell['x'] = 0		# Координата икс.
	dbell['y'] = 0 		
	dbell['z'] = 0
	
	dbell['glob_g'] = 9.86
	
	
	# Стартовый участок.
	print("WNDTWR = Стартовый участок.")
	import let_point_thy_path as willowisp
	# Если ракета может навестись на цель сразу, то сразу же и наводим.
	if(ddata['if_povorot_start'] == 1):
		#Здесь вычисляем угол поворота. 
		print("WNDTWR = Вычисление угла поворота.")
		
		teta_gorizont = willowisp.water_teta(xc, zc, 0, 0)
		teta_vertical = willowisp.water_teta(xc, yc, 0, 0)
		dbell['teta_sumxy'] = teta_vertical
		dbell['teta_sumxz'] = teta_gorizont
		print("\tПоворот выполнен.")
		print('\t\t', dbell["teta_sumxy"])
		print('\t\t', dbell['teta_sumxz'])
	
	print("WNDTWR = Расчёт стартового участка.")
	while(m_0fuel > 0):
		# 1. Получить текущую плотность воздушных масс.
		rho = atm.f_rho_dra(y, dra)
		dbell['rho'] = rho
		# 2. Провести расчёт.
		t += dt
		m_0fuel -= G_10*dt
		dbell = bellkeeper(dbell, ldc, rho, dt=ddata['dt'], flag_G=False)
		# 3. Снятие контрольных данных в списки.
		l_m.append(dbell['m'])
		#l_vxy.append(dbell['V_xy'])
		#l_vyy.append(dbell['V_yy'])
		#l_vzx.append(dbell['V_zx'])
		#l_vzz.append(dbell['V_zz'])
		l_vsum.append(dbell['V_sum'])
		l_t.append(t)
		l_x.append(dbell['x'])
		l_y.append(dbell['y'])
		l_z.append(dbell['z'])
		# 4. Блок проверки флагов.
		L = L_to_point(dbell['x'], dbell['y'], dbell['z'], xc, yc, zc)
		#print("L = ", L)
		#print("xc = %s, yc = %s, zc = %s" %(xc, yc, zc))
		if(L <= epsilon):
			# Попадание.
			answer = inf_answer("start")
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			# Здесь нужно добавить рисование графиков в точку.
			return(dreturn)
		
		lr={}
		lr['x'] = dbell['x']
		lr['y'] = dbell['y']
		lr['z'] = dbell['z']
		#lr['V_xy'] = dbell['V_xy']
		#lr['V_yy'] = dbell['V_yy']
		#lr['V_zx'] = dbell['V_zx']
		#lr['V_zz'] = dbell['V_zz']
		lr['V_sum'] = dbell['V_sum']
		#lr['psi'] = dbell['teta_sumxz']
		#lr['nu'] = dbell['teta_sumxy']
		lr['psi'] = dbell['psi']
		lr['nu'] = dbell['nu_v']
		lr['gamma'] = dbell['gamma']
		lr['limpsi'] = limpsi
		lr['limnu'] = limnu
		lr['limgamma'] = limgamma
		
		lc={}
		lc['x'] = xc
		lc['y'] = yc
		lc['z'] = zc
		lc['V_xy'] = 0
		lc['V_yy'] = 0
		lc['V_xz'] = 0
		lc['V_zz'] = 0
		lc['psi'] = 0
		lc['nu'] = 0
		letter = 'to point'
		if(dbell['V_sum'] > ddata['V_manevr'] and ddata['flag_startmanevr'] != 0):
			# Вычисление угла поворота по методу наведения.
			mg = dbell['m'] * dbell['glob_g']
			dget = willowisp.guilding_star(dbell['P'], dbell['X'], dbell['Y'], mg, lr, lc, letter, flag_recommend=True, teta_recommend=2.5)
			# Проверка.
			if(dget == False):
				answer = inf_answer("teta")
				dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
				return(dreturn)
			#dbell['teta_xz'] = dget['teta_xz']
			#dbell['teta_xy'] = dget['teta_xy']
			#dbell['teta_sumxz'] = dget['teta_sumxz']
			#dbell['teta_sumxy'] = dget['teta_sumxy']
			#dbell['gamma'] = dget['gamma']
			dbell['nu_p'] = dget['nu'] + 2.1
			dbell['nu_v'] = dget['nu']
			dbell['psi'] = dget['psi']
			
			# Пифагоровы штаны во все стороны равны.
			dbell['V_x'] = dbell['V_sum'] * cos(dbell['nu_v'])
			dbell['V_y'] = dbell['V_sum'] * sin(dbell['nu_v'])
			
			# Подмена угла для реализации нужного полёта.
			#dbell['teta_xz'] = 0
			#dbell['teta_xy'] = 3
			#dbell['teta_sumxz'] = teta_gorizont
			#dbell['teta_sumxy'] = teta_vertical
			
			flag = False
			for elem in dbell:
				if(elem == None):
					flag = True
			if(flag == True):
				print("Флажок поднят, углы не определены.")
			
		# Поправка атмосферы.
		y = lr['y']
		if(y > altmax):
			answer = inf_answer("atmosphere")
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			return(dreturn)
	
	# Переходный процесс. Старт -> Марш.
	print("WNDTWR = переходный процесс.")
	m_sbros = ddata['m_0'] - ddata['m_0fuel']
	print("m_sbros = ", m_sbros)
	if(ddata['if_sbros_march'] == 1):
		dbell['m'] -= m_sbros
		l_m.append(dbell['m'])
		#l_vxy.append(dbell['V_xy'])
		#l_vyy.append(dbell['V_yy'])
		#l_vzx.append(dbell['V_zx'])
		#l_vzz.append(dbell['V_zz'])
		l_vsum.append(dbell['V_sum'])
		l_t.append(t)
		l_x.append(dbell['x'])
		l_y.append(dbell['y'])
		l_z.append(dbell['z'])
	# Проверка по набору скорости в конце стартового участка.
	M = dbell['V_sum']/360
	if(M < ddata['M_start_end']):
		# Не набрана достаточная скорость..
		answer = inf_answer("M<V")
		dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
		# Здесь нужно добавить рисование графиков в точку.
		return(dreturn)
	
	
	# Маршевый участок.
	print("WNDTWR = маршевый участок.")
	# Проверки поворота устраивать не буду.
	while(m_1fuel > 0):
		#print("\tm_1fuel = ", m_1fuel)
		#print("\tV_x = ", l_vx[-1])
		#print("\tV_z = ", l_vz[-1])
		#print("\tV_y = ", l_vy[-1])
		print("V_sum = ", dbell['V_sum'])
		#if(m_1fuel <= 560):
		#	# Попадание.
		#	answer = inf_answer("start")
		#	dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vx, l_vy, l_vz, l_m, l_resist, answer, True)
		#	# Здесь нужно добавить рисование графиков в точку.
		#	return(dreturn)
		# 1. Получить текущую плотность воздушных масс.
		rho = atm.f_rho_dra(y, dra)
		dbell['rho'] = rho
		# 2. Провести расчёт.
		t += dt
		m_1fuel -= G_11*dt
		dbell = bellkeeper(dbell, ldc, rho, dt=ddata['dt'], flag_G=True, flag_P=True)
		# 3. Снятие контрольных данных в списки.
		l_m.append(dbell['m'])
		#l_vxy.append(dbell['V_xy'])
		#l_vyy.append(dbell['V_yy'])
		#l_vzx.append(dbell['V_zx'])
		#l_vzz.append(dbell['V_zz'])
		l_vsum.append(dbell['V_sum'])
		l_t.append(t)
		l_x.append(dbell['x'])
		l_y.append(dbell['y'])
		l_z.append(dbell['z'])
		# 4. Блок проверки флагов.
		L = L_to_point(dbell['x'], dbell['y'], dbell['z'], xc, yc, zc)
		print("L = ", L)
		#print("xc = %s, yc = %s, zc = %s" %(xc, yc, zc))
		if(L <= epsilon):
			# Попадание.
			answer = inf_answer("march")
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			# Здесь нужно добавить рисование графиков в точку.
			return(dreturn)
		#else:
		#	# Промах.
		#	answer = inf_answer("fall")
		#	dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vx, l_vy, l_vz, l_m, l_resist, answer, False)
		#	# Здесь нужно добавить рисование графиков в точку.
		#	return(dreturn)
		# 5. Блок реализации поворота.
		lr={}
		lr['x'] = dbell['x']
		lr['y'] = dbell['y']
		lr['z'] = dbell['z']
		lr['V_x'] = dbell['V_x']
		lr['V_y'] = dbell['V_y']
		#lr['V_xy'] = dbell['V_xy']
		#lr['V_yy'] = dbell['V_yy']
		#lr['V_zx'] = dbell['V_zx']
		#lr['V_zz'] = dbell['V_zz']
		#lr['psi'] = dbell['teta_sumxz']
		#lr['nu'] = dbell['teta_sumxy']
		lr['psi'] = dbell['psi']
		lr['nu'] = dbell['nu_v']
		lr['gamma'] = dbell['gamma']
		lr['limpsi'] = limpsi
		lr['limnu'] = limnu
		lr['limgamma'] = limgamma
		
		lc={}
		lc['x'] = xc
		lc['y'] = yc
		lc['z'] = zc
		lc['V_x'] = 0
		lc['V_y'] = 0
		lc['V_z'] = 0
		lc['psi'] = 0
		lc['nu'] = 0
		letter = 'to point'
		
		# Вычисление угла поворота по методу наведения.
		mg = dbell['m'] * dbell['glob_g']
		dget = willowisp.guilding_star(dbell['P'], dbell['X'], dbell['Y'], mg, lr, lc, letter, flag_recommend=True, teta_recommend=2.5)
		# Проверка.
		if(dget == False):
			answer = inf_answer("teta")
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			return(dreturn)
		
		#dbell['teta_xz'] = dget['teta_xz']
		#dbell['teta_xy'] = dget['teta_xy']
		#dbell['teta_sumxz'] = dget['teta_sumxz']
		#dbell['teta_sumxy'] = dget['teta_sumxy']
		#dbell['gamma'] = dget['gamma']
		
		dbell['nu_v'] = dget['nu'] 
		dbell['nu_p'] = dget['nu'] + 2.2
		dbell['psi'] = dget['psi']
		
		# Пифагоровы штаны во все стороны равны.
		dbell['V_x'] = dbell['V_sum'] * cos(dbell['nu_v'])
		dbell['V_y'] = dbell['V_sum'] * sin(dbell['nu_v'])
		
		# Поправка атмосферы.
		y = lr['y']
		if(y > altmax):
			answer = inf_answer("atmosphere")
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			return(dreturn)
		del lr, lc
	
	
	# Пассивный полёт. И никаких вам переходных процессов.
	print("WNDTWR = пассивный полёт.")
	answer = inf_answer("teta")
	dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
	return(dreturn)
	dbell['P_xy'] = 0
	dbell['P_xz'] = 0
	dbell['G_10'] = 0
	dbell['G_11'] = 0
	
	V = dbell['V_sum']
	dbell['P_10'] = 0
	dbell['P_11'] = 0
	counter = 0
	while(V > ddata['M_passive_end'] and dbell['nu_p'] > -90):
		#print("\tV = ", V)
		#print('\tM_passive_end = ', ddata['M_passive_end'])
		# 1. Получить текущую плотность воздушных масс.
		rho = atm.f_rho_dra(y, dra)
		dbell['rho'] = rho
		# 2. Провести расчёт.
		t += dt
		dbell = bellkeeper(dbell, ldc, rho, dt=ddata['dt'])
		# 3. Снятие контрольных данных в списки.
		l_m.append(dbell['m'])
		#l_vxy.append(dbell['V_xy'])
		#l_vyy.append(dbell['V_yy'])
		#l_vzx.append(dbell['V_zx'])
		#l_vzz.append(dbell['V_zz'])
		l_vsum.append(dbell['V_sum'])
		l_t.append(t)
		l_x.append(dbell['x'])
		l_y.append(dbell['y'])
		l_z.append(dbell['z'])
		# 4. Блок проверки флагов.
		L = L_to_point(dbell['x'], dbell['y'], dbell['z'], x, y, z)
		print("L = ", L)
		if(L <= epsilon):
			# Попадание.
			answer = inf_answer("fatigue")
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			# Здесь нужно добавить рисование графиков в точку.
			return(dreturn)
		#else:
		#	# Промах.
		#	answer = inf_answer("M<V")
		#	dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vx, l_vy, l_vz, l_m, l_resist, answer, False)
		#	# Здесь нужно добавить рисование графиков в точку.
		#	return(dreturn)
		# 5. Блок реализации поворота.
		lr={}
		lr['x'] = dbell['x']
		lr['y'] = dbell['y']
		lr['z'] = dbell['z']
		lr['V_x'] = dbell['V_x']
		lr['V_y'] = dbell['V_y']
		#lr['V_xy'] = dbell['V_xy']
		#lr['V_yy'] = dbell['V_yy']
		#lr['V_zx'] = dbell['V_zx']
		#lr['V_zz'] = dbell['V_zz']
		lr['V_sum'] = dbell['V_sum']
		#lr['psi'] = dbell['teta_sumxz']
		#lr['nu'] = dbell['teta_sumxy']
		lr['psi'] = dbell['psi']
		lr['nu'] = dbell['nu_v']
		lr['gamma'] = dbell['gamma']
		lr['limpsi'] = limpsi
		lr['limnu'] = limnu
		lr['limgamma'] = limgamma
		lc={}
		lc['x'] = xc
		lc['y'] = yc
		lc['z'] = zc
		lc['V_x'] = 0
		lc['V_y'] = 0
		lc['V_z'] = 0
		lc['psi'] = 0
		
		# Вычисление угла поворота по методу наведения.
		mg = dbell['m'] * dbell['glob_g']
		dget = willowisp.guilding_star(0, dbell['X'], dbell['Y'], mg, lr, lc, letter)
		# Проверка.
		if(dget == False):
			print("Пассивный полёт, ракета уходит в сторону?")
			answer = inf_answer("teta")
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			return(dreturn)
		#dbell['teta_xz'] = dget['teta_xz']
		#dbell['teta_xy'] = dget['teta_xy']
		#dbell['teta_sumxz'] = dget['teta_sumxz']
		#dbell['teta_sumxy'] = dget['teta_sumxy']
		#dbell['gamma'] = dget['gamma']
		dbell['nu_v'] = dget['nu']
		dbell['nu_p'] = dget['nu'] + 2
		dbell['psi'] = dget['psi']
		
		V = dbell['V_sum']
		M = V/360
		
		# Поправка атмосферы.
		y = lr['y']
		if(y > altmax):
			answer = inf_answer("atmosphere")
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			return(dreturn)
		del lr, lc
		
	
	# Наше присутсвие здесь значит собой промах ракеты. Не знаю, почему.
	answer = inf_answer("fall")
	dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
	return(dreturn)











































































#	Общие расчёты по пространственной задаче при полёте в точку.
#	Получает:			1) ddata - словарь данных,
#						2) x, y, z - координаты цели.
# 						3) flag_graph - Флаг формы возрата решения.
# 						4) flag_ret_all - Флаг возврата данных решения.
#	Делает:				Интегрирует.
#	Возвращает:			1)	Если флаг <flag_graph> опущен:
# 								Печать в консоль, сбор данных в контейнер и его возврат.
# 							Если поднят:
# 								Строит графики.
# 						2)	Если флаг <flag_ret_all> опущен:
# 								Возврат <False> если решение неудовлетворяет.
# 								Возврат <True> если решение удовлетворительно.
# 							Если флаг <flag_ret_all> поднят:
# 								Возврат словаря из вычисленных списков. 
#	Описание:	Процедура должна получить уже сформированный словарь характеристик ddata. Процедура возвращает
#	обратно не словарь <ddata>, а словарь <dret>, который содержит:
# 										- l_x = [] - список координат горизонтали.
# 										- l_y = [] - список координат вертикали.
# 										- l_t = [] - список координат времени.
# 										- l_V = [] - список скоростей.
# 										- l_m = [] - список масс.
# 										- l_resist = [] - список сопротивлений.
# 										- answer - строка состояния решения. Доступны следующие состояния:
# 												"Ракета достигла цели на стартовом участке."
# 												"Ракета достигла цели на маршевом участке."
# 												"Ракета достигла цели на пассивном участке."
# 												"Ракета не долетела."
def d_to_ALL(letter, ddata, xc, yc, zc, flag_graph=False, flag_ret_all=False):
	from trigonometria import sin, cos
	import atmosphere_GOST4401_81 as atm
	import cxyz 
	def inf_dreturn(l_x, l_y, l_z, l_t, l_V_xy, l_V_yy, l_V_zx, l_V_zz, l_V_sum, l_m, answer, flag):
		dreturn = {}
		dreturn['l_x'] = l_x.copy()
		dreturn['l_y'] = l_y.copy()
		dreturn['l_z'] = l_z.copy()
		dreturn['l_t'] = l_t.copy()
		dreturn['l_V_xy'] = l_V_xy.copy()
		dreturn['l_V_yy'] = l_V_yy.copy()
		dreturn['l_V_zx'] = l_V_zx.copy()
		dreturn['l_V_zz'] = l_V_zz.copy()
		dreturn['l_V_sum'] = l_V_sum.copy()
		dreturn['l_m'] = l_m.copy()
		dreturn['answer'] = answer
		dreturn['flag'] = flag
		return(dreturn)
	
	# Формовка строки характеризующей решение:
	def inf_answer(position, retflag=False):
		if(position=='start'):
			if(retflag):
				return(True)
			return("Ракета достигла своей цели на стартовом участке.")
		if(position=='march'):
			if(retflag):
				return(True)
			return('Ракета достигла своей цели на маршевом участке.')
		if(position=='fatigue'):
			if(retflag):
				return(True)
			return('Ракета достигла своей цели на пассивном участке.')
		if(position=='fall'):
			if(retflag):
				return(False)
			return('Ракета не долетела.')
		if(position=="M<V"):
			if(retflag):
				return(False)
			return('Ракета не набрала нужной скорости.')
		if(position=="teta"):
			if(retflag):
				return(False)
			return("Ракета начала вести угол атаки более 10 градусов, что не есть хорошо.")
		if(position=='atmosphere'):
			if(retflag):
				return(False)
			return("Ракета вылетела за пределы атмосферы.")
		else:
			if(retflag):
				return(False)
			return('Здесь был бред')
		
	# Функция расчёта дальности до цели.
	def L_to_point(x, y, z, xp, yp, zp):
		xt = abs(xp - x)
		yt = abs(yp - y)
		zt = abs(zp - z)
		#print('\txt = ',xt)
		#print('\tyt = ',yt)
		#print('\tzt = ',zt)
		L = (xt**2 + yt**2 + zt**2)**0.5
		return(L)
	
    	# Перед решением получаем начальные данные.
	print("WNDTWR= получение начальных данных.")
	#dra = rho()
	dra = atm.get_table(fcsv=True)	#	Таблица плотностей атмосферы в зависимости от высоты.
	altmax = atm.getmax()	#	Получили таблицу зависимости атмосферы.
	dcx = {}	# 	Получили таблицу зависимости коэффициентов Cx, Cy, Cz.
	dcy = {}
	dcz = {}
	dcx['lx'], dcx['ly'] = cxyz.get_table("cx")
	dcy['lx'], dcy['ly'] = cxyz.get_table("cy")
	dcz['lx'], dcz['ly'] = cxyz.get_table("cz")
	ldc = []
	ldc.append(dcx)
	ldc.append(dcy)
	ldc.append(dcz)
	
	
    	
	# Перед решением снимаем копию данных для удобства обращения.
	#print(ddata)
	mu_0 = ddata['mu_0']
	mu_1 = ddata['mu_1']
	d = ddata['d']
	F_midelya = ddata['F_midelya']
	m_sum = ddata['m_sum']
	m_0 = ddata['m_0']
	m_0fuel = ddata['m_0fuel']
	m_1 = ddata['m_1']
	m_1fuel = ddata['m_1fuel']
	P_10 = ddata['P_10']
	P_11 = ddata['P_11']
	G_10 = ddata['G_10']
	G_11 = ddata['G_11']
	limpsi = ddata['limpsi']
	limnu = ddata['limnu']
	limgamma = ddata['limgamma']
	
	# Принимаем начальные данные.
	dt = ddata['dt']		# Шаг интегрирования по времени.
	t = 0			# Общее время.
	epsilon_m = 5	# Ошибка определения масс топлива.
	x = 0			# Начальная горизонтальная координата <x>.
	y = 0			# Начальная вертикальная координата.
	z = 0			# Начальная координата горизонта <z>.
	V = 0			# Начальная скорость.
	epsilon = 30	# Вероятное отклонение от цели, засчитываемое за попадание.
    	
	l_m = []
	l_vxy = []
	l_vyy = []
	l_vzx = []
	l_vzz = []
	l_vsum = []
	l_t = []
	l_x = []
	l_y = []
	l_z = []
	l_m.append(m_sum)
	l_vxy.append(0)
	l_vyy.append(0)
	l_vzx.append(0)
	l_vzz.append(0)
	l_vsum.append(0)
	l_t.append(0)
	l_x.append(0)
	l_y.append(0)
	l_z.append(0)
	
	# Теперь нужно сформировать словарь, который будет гоняться в функции <bellkeeper>.
	dbell = {}
	
	# Скорости.
	#dbell['V_xy'] = 0	# Скорости поступательного движения по осям.
	#dbell['V_yy'] = 0
	#dbell['V_zx'] = 0
	#dbell['V_zz'] = 0
	dbell['V_x'] = 0
	dbell['V_y'] = 0
	dbell['V_sum'] = 0
	
	# Массы.
	dbell['m'] = m_sum	# Масса
	dbell['m_fuel'] = m_0fuel + m_1fuel
	
	# Тяги.
	dbell['P'] = P_10	# Значения тяги.
	dbell['P_10'] = P_10
	dbell['P_11'] = P_11
	dbell['G_10'] = G_10
	dbell['G_11'] = G_11
	
	# Углы полёта.
	dbell['teta_sumxy'] = ddata['teta_0']	# Суммарный угол старта ракеты в вертикальной плоскости.
	dbell['teta_sumxz'] = ddata['psi']		# Суммарный угол старта ракеты в горизонтальной плоскости.
	dbell['teta_xy'] = 3	# Значение угла скольжения ЛА. Вертикальная плоскость.
	dbell['teta_xz'] = 0	# Значения угла скольжения ЛА в горизонтальной плоскости.
	dbell['omega_x'] = 0	# Угловая скорость.
	dbell['omega_y'] = 0	
	dbell['omega_z'] = 0
	dbell['gamma'] = ddata['gamma']		# Угол крена стартовый, обычно, равен нулю.
	
	dbell['nu_p'] = ddata['teta_0']
	dbell['nu_v'] = ddata['teta_0']
	dbell['psi'] = ddata['psi']
	
	# Значения АД-сил.
	dbell['X_xy'] = 0		# Значение АД-сил.
	dbell['Y_xy'] = 0
	dbell['X_xz'] = 0
	dbell['Z_xz'] = 0
	dbell['X'] = 0
	dbell['Y'] = 0
	
	# Значения моментов инерции и прочего.
	dbell['I_x'] = ddata['I_x']	# Моменты инерции вдоль осей движения.
	dbell['I_y'] = ddata['I_y']
	dbell['I_z'] = ddata['I_z']
	dbell['M_x'] = ddata['M_x']	# Балансировочные зависимости ЛА.
	dbell['M_y'] = ddata['M_y']
	dbell['M_z'] = ddata['M_z']
	
	dbell['x'] = 0		# Координата икс.
	dbell['y'] = 0 		
	dbell['z'] = 0
	
	dbell['glob_g'] = 9.86
	
	# Цель.
	vcx = ddata['v_0x']
	vcy = ddata['v_0y']
	vcz = ddata['v_0z']
	
	
	# Стартовый участок.
	#print("WNDTWR = Стартовый участок.")
	import let_point_thy_path as willowisp
	# Если ракета может навестись на цель сразу, то сразу же и наводим.
	if(ddata['if_povorot_start'] == 1):
		#Здесь вычисляем угол поворота. 
		#print("WNDTWR = Вычисление угла поворота.")
		
		teta_gorizont = willowisp.water_teta(xc, zc, 0, 0)
		teta_vertical = willowisp.water_teta(xc, yc, 0, 0)
		dbell['teta_sumxy'] = teta_vertical
		dbell['teta_sumxz'] = teta_gorizont
		#print("\tПоворот выполнен.")
		#print('\t\t', dbell["teta_sumxy"])
		#print('\t\t', dbell['teta_sumxz'])
	
	#print("WNDTWR = Расчёт стартового участка.")
	while(m_0fuel > 0):
		# 1. Получить текущую плотность воздушных масс.
		rho = atm.f_rho_dra(y, dra)
		dbell['rho'] = rho
		# 2. Провести расчёт.
		t += dt
		m_0fuel -= G_10*dt
		dbell = bellkeeper(dbell, ldc, rho, dt=ddata['dt'], flag_G=False)
		# 3. Снятие контрольных данных в списки.
		l_m.append(dbell['m'])
		#l_vxy.append(dbell['V_xy'])
		#l_vyy.append(dbell['V_yy'])
		#l_vzx.append(dbell['V_zx'])
		#l_vzz.append(dbell['V_zz'])
		l_vsum.append(dbell['V_sum'])
		l_t.append(t)
		l_x.append(dbell['x'])
		l_y.append(dbell['y'])
		l_z.append(dbell['z'])
		# 4. Блок проверки флагов.
		L = L_to_point(dbell['x'], dbell['y'], dbell['z'], xc, yc, zc)
		#print("L = ", L)
		#print("xc = %s, yc = %s, zc = %s" %(xc, yc, zc))
		if(L <= epsilon):
			# Попадание.
			answer = inf_answer("start", retflag=True)
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			# Здесь нужно добавить рисование графиков в точку.
			return(dreturn)
		
		lr={}
		lr['x'] = dbell['x']
		lr['y'] = dbell['y']
		lr['z'] = dbell['z']
		#lr['V_xy'] = dbell['V_xy']
		#lr['V_yy'] = dbell['V_yy']
		#lr['V_zx'] = dbell['V_zx']
		#lr['V_zz'] = dbell['V_zz']
		lr['V_sum'] = dbell['V_sum']
		#lr['psi'] = dbell['teta_sumxz']
		#lr['nu'] = dbell['teta_sumxy']
		lr['psi'] = dbell['psi']
		lr['nu'] = dbell['nu_v']
		lr['gamma'] = dbell['gamma']
		lr['limpsi'] = limpsi
		lr['limnu'] = limnu
		lr['limgamma'] = limgamma
		
		lc={}
		lc['x'] = xc
		lc['y'] = yc
		lc['z'] = zc
		lc['V_xy'] = 0
		lc['V_yy'] = 0
		lc['V_xz'] = 0
		lc['V_zz'] = 0
		lc['psi'] = 0
		lc['nu'] = 0
		if(dbell['V_sum'] > ddata['V_manevr'] and ddata['flag_startmanevr'] != 0):
			# Вычисление угла поворота по методу наведения.
			mg = dbell['m'] * dbell['glob_g']
			dget = willowisp.guilding_star(dbell['P'], dbell['X'], dbell['Y'], mg, lr, lc, letter, flag_recommend=True, teta_recommend=2.5)
			# Проверка.
			if(dget == False):
				answer = inf_answer("teta")
				dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
				return(dreturn)
			#dbell['teta_xz'] = dget['teta_xz']
			#dbell['teta_xy'] = dget['teta_xy']
			#dbell['teta_sumxz'] = dget['teta_sumxz']
			#dbell['teta_sumxy'] = dget['teta_sumxy']
			#dbell['gamma'] = dget['gamma']
			dbell['nu_p'] = dget['nu'] + 2.1
			dbell['nu_v'] = dget['nu']
			dbell['psi'] = dget['psi']
			
			# Пифагоровы штаны во все стороны равны.
			dbell['V_x'] = dbell['V_sum'] * cos(dbell['nu_v'])
			dbell['V_y'] = dbell['V_sum'] * sin(dbell['nu_v'])
			
			# Подмена угла для реализации нужного полёта.
			#dbell['teta_xz'] = 0
			#dbell['teta_xy'] = 3
			#dbell['teta_sumxz'] = teta_gorizont
			#dbell['teta_sumxy'] = teta_vertical
			
			flag = False
			for elem in dbell:
				if(elem == None):
					flag = True
			if(flag == True):
				print("Флажок поднят, углы не определены.")
			
		# Поправка атмосферы.
		y = lr['y']
		if(y > altmax-10000 or y <0):
			answer = inf_answer("atmosphere", retflag=True)
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			return(dreturn)
		
		# Полёт цели.
		xc += vcx*dt
		yc += vcy*dt
		zc += vcz*dt
	
	# Переходный процесс. Старт -> Марш.
	#print("WNDTWR = переходный процесс.")
	m_sbros = ddata['m_0'] - ddata['m_0fuel']
	#print("m_sbros = ", m_sbros)
	if(ddata['if_sbros_march'] == 1):
		dbell['m'] -= m_sbros
		l_m.append(dbell['m'])
		#l_vxy.append(dbell['V_xy'])
		#l_vyy.append(dbell['V_yy'])
		#l_vzx.append(dbell['V_zx'])
		#l_vzz.append(dbell['V_zz'])
		l_vsum.append(dbell['V_sum'])
		l_t.append(t)
		l_x.append(dbell['x'])
		l_y.append(dbell['y'])
		l_z.append(dbell['z'])
	# Проверка по набору скорости в конце стартового участка.
	M = dbell['V_sum']/360
	if(M < ddata['M_start_end']):
		# Не набрана достаточная скорость..
		answer = inf_answer("M<V", retflag=True)
		dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
		# Здесь нужно добавить рисование графиков в точку.
		return(dreturn)
	
	
	# Маршевый участок.
	#print("WNDTWR = маршевый участок.")
	# Проверки поворота устраивать не буду.
	while(m_1fuel > 0):
		#print("\tm_1fuel = ", m_1fuel)
		#print("\tV_x = ", l_vx[-1])
		#print("\tV_z = ", l_vz[-1])
		#print("\tV_y = ", l_vy[-1])
		#print("V_sum = ", dbell['V_sum'])
		#if(m_1fuel <= 560):
		#	# Попадание.
		#	answer = inf_answer("start")
		#	dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vx, l_vy, l_vz, l_m, l_resist, answer, True)
		#	# Здесь нужно добавить рисование графиков в точку.
		#	return(dreturn)
		# 1. Получить текущую плотность воздушных масс.
		rho = atm.f_rho_dra(y, dra)
		dbell['rho'] = rho
		# 2. Провести расчёт.
		t += dt
		m_1fuel -= G_11*dt
		dbell = bellkeeper(dbell, ldc, rho, dt=ddata['dt'], flag_G=True, flag_P=True)
		# 3. Снятие контрольных данных в списки.
		l_m.append(dbell['m'])
		#l_vxy.append(dbell['V_xy'])
		#l_vyy.append(dbell['V_yy'])
		#l_vzx.append(dbell['V_zx'])
		#l_vzz.append(dbell['V_zz'])
		l_vsum.append(dbell['V_sum'])
		l_t.append(t)
		l_x.append(dbell['x'])
		l_y.append(dbell['y'])
		l_z.append(dbell['z'])
		# 4. Блок проверки флагов.
		L = L_to_point(dbell['x'], dbell['y'], dbell['z'], xc, yc, zc)
		#print("L = ", L)
		#print("xc = %s, yc = %s, zc = %s" %(xc, yc, zc))
		if(L <= epsilon):
			# Попадание.
			answer = inf_answer("march", retflag=True)
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			# Здесь нужно добавить рисование графиков в точку.
			return(dreturn)
		#else:
		#	# Промах.
		#	answer = inf_answer("fall")
		#	dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vx, l_vy, l_vz, l_m, l_resist, answer, False)
		#	# Здесь нужно добавить рисование графиков в точку.
		#	return(dreturn)
		# 5. Блок реализации поворота.
		lr={}
		lr['x'] = dbell['x']
		lr['y'] = dbell['y']
		lr['z'] = dbell['z']
		lr['V_x'] = dbell['V_x']
		lr['V_y'] = dbell['V_y']
		#lr['V_xy'] = dbell['V_xy']
		#lr['V_yy'] = dbell['V_yy']
		#lr['V_zx'] = dbell['V_zx']
		#lr['V_zz'] = dbell['V_zz']
		#lr['psi'] = dbell['teta_sumxz']
		#lr['nu'] = dbell['teta_sumxy']
		lr['psi'] = dbell['psi']
		lr['nu'] = dbell['nu_v']
		lr['gamma'] = dbell['gamma']
		lr['limpsi'] = limpsi
		lr['limnu'] = limnu
		lr['limgamma'] = limgamma
		
		lc={}
		lc['x'] = xc
		lc['y'] = yc
		lc['z'] = zc
		lc['V_x'] = 0
		lc['V_y'] = 0
		lc['V_z'] = 0
		lc['psi'] = 0
		lc['nu'] = 0
		letter = 'to point'
		
		# Вычисление угла поворота по методу наведения.
		mg = dbell['m'] * dbell['glob_g']
		dget = willowisp.guilding_star(dbell['P'], dbell['X'], dbell['Y'], mg, lr, lc, letter, flag_recommend=True, teta_recommend=2.5)
		# Проверка.
		if(dget == False):
			answer = inf_answer("teta", retflag=True)
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			return(dreturn)
		
		#dbell['teta_xz'] = dget['teta_xz']
		#dbell['teta_xy'] = dget['teta_xy']
		#dbell['teta_sumxz'] = dget['teta_sumxz']
		#dbell['teta_sumxy'] = dget['teta_sumxy']
		#dbell['gamma'] = dget['gamma']
		
		dbell['nu_v'] = dget['nu'] 
		dbell['nu_p'] = dget['nu'] + 2.2
		dbell['psi'] = dget['psi']
		
		# Пифагоровы штаны во все стороны равны.
		dbell['V_x'] = dbell['V_sum'] * cos(dbell['nu_v'])
		dbell['V_y'] = dbell['V_sum'] * sin(dbell['nu_v'])
		
		# Поправка атмосферы.
		y = lr['y']
		if(y > altmax-10000 or y <0):
			answer = inf_answer("atmosphere", retflag=True)
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			return(dreturn)
		del lr, lc
		
		
		# Полёт цели.
		xc += vcx*dt
		yc += vcy*dt
		zc += vcz*dt
		#print("\t\tЦель:\n\t\txc = %s\n\t\tyc= %s\n\t\tzc = %s" %(xc, yc, zc))
	
	
	# Пассивный полёт. И никаких вам переходных процессов.
	#print("WNDTWR = пассивный полёт.")
	#answer = inf_answer("teta")
	#dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
	#return(dreturn)
	dbell['P_xy'] = 0
	dbell['P_xz'] = 0
	dbell['G_10'] = 0
	dbell['G_11'] = 0
	
	V = dbell['V_sum']
	dbell['P_10'] = 0
	dbell['P_11'] = 0
	counter = 0
	while(V > ddata['M_passive_end'] and dbell['nu_p'] > -90):
		#print("\tV = ", V)
		#print('\tM_passive_end = ', ddata['M_passive_end'])
		# 1. Получить текущую плотность воздушных масс.
		rho = atm.f_rho_dra(y, dra)
		dbell['rho'] = rho
		# 2. Провести расчёт.
		t += dt
		dbell = bellkeeper(dbell, ldc, rho, dt=ddata['dt'])
		# 3. Снятие контрольных данных в списки.
		l_m.append(dbell['m'])
		#l_vxy.append(dbell['V_xy'])
		#l_vyy.append(dbell['V_yy'])
		#l_vzx.append(dbell['V_zx'])
		#l_vzz.append(dbell['V_zz'])
		l_vsum.append(dbell['V_sum'])
		l_t.append(t)
		l_x.append(dbell['x'])
		l_y.append(dbell['y'])
		l_z.append(dbell['z'])
		# 4. Блок проверки флагов.
		L = L_to_point(dbell['x'], dbell['y'], dbell['z'], xc, yc, zc)
		#print("L = ", L)
		if(L <= epsilon):
			# Попадание.
			answer = inf_answer("fatigue", retflag=True)
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			# Здесь нужно добавить рисование графиков в точку.
			return(dreturn)
		#else:
		#	# Промах.
		#	answer = inf_answer("M<V")
		#	dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vx, l_vy, l_vz, l_m, l_resist, answer, False)
		#	# Здесь нужно добавить рисование графиков в точку.
		#	return(dreturn)
		# 5. Блок реализации поворота.
		lr={}
		lr['x'] = dbell['x']
		lr['y'] = dbell['y']
		lr['z'] = dbell['z']
		lr['V_x'] = dbell['V_x']
		lr['V_y'] = dbell['V_y']
		#lr['V_xy'] = dbell['V_xy']
		#lr['V_yy'] = dbell['V_yy']
		#lr['V_zx'] = dbell['V_zx']
		#lr['V_zz'] = dbell['V_zz']
		lr['V_sum'] = dbell['V_sum']
		#lr['psi'] = dbell['teta_sumxz']
		#lr['nu'] = dbell['teta_sumxy']
		lr['psi'] = dbell['psi']
		lr['nu'] = dbell['nu_v']
		lr['gamma'] = dbell['gamma']
		lr['limpsi'] = limpsi
		lr['limnu'] = limnu
		lr['limgamma'] = limgamma
		lc={}
		lc['x'] = xc
		lc['y'] = yc
		lc['z'] = zc
		lc['V_x'] = 0
		lc['V_y'] = 0
		lc['V_z'] = 0
		lc['psi'] = 0
		
		# Вычисление угла поворота по методу наведения.
		mg = dbell['m'] * dbell['glob_g']
		dget = willowisp.guilding_star(0, dbell['X'], dbell['Y'], mg, lr, lc, letter)
		# Проверка.
		if(dget == False):
			print("Пассивный полёт, ракета уходит в сторону?")
			answer = inf_answer("teta", retflag=True)
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			return(dreturn)
		#dbell['teta_xz'] = dget['teta_xz']
		#dbell['teta_xy'] = dget['teta_xy']
		#dbell['teta_sumxz'] = dget['teta_sumxz']
		#dbell['teta_sumxy'] = dget['teta_sumxy']
		#dbell['gamma'] = dget['gamma']
		dbell['nu_v'] = dget['nu']
		dbell['nu_p'] = dget['nu'] + 2
		dbell['psi'] = dget['psi']
		
		V = dbell['V_sum']
		M = V/360
		
		# Поправка атмосферы.
		y = lr['y']
		if(y > altmax-10000 or y < 0):
			answer = inf_answer("atmosphere", retflag=True)
			dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
			return(dreturn)
		del lr, lc
		
		
		# Полёт цели.
		xc += vcx*dt
		yc += vcy*dt
		zc += vcz*dt
		
	
	# Наше присутсвие здесь значит собой промах ракеты. Не знаю, почему.
	answer = inf_answer("fall", retflag=True)
	dreturn = inf_dreturn(l_x, l_y, l_z, l_t, l_vxy, l_vyy, l_vzx, l_vzz, l_vsum, l_m, answer, True)
	return(dreturn)
