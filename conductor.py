#	
#	
#	
#	
#	
#	ЗДЕСЬ ТОЛЬКО НАВЕДЕНИЕ И БОЛЕЕ НИЧЕГО.
#	
#	
#	
#	
#	
import math


#	=	Объёмный мир	=

# Процедура вычисляет угол линии по двум точкам к горизонту. Для плоской задачи.
# Получает:                                     4 координаты:
#                                                               -       xc              - координата x цели.
#                                                               -       yc              - координата y цели.
#                                                               -       xr              - координата x ракеты.
#                                                               -       yr              - координата y ракеты.
# Делает:                                       Вычисляет угол.
# Возвращает:                           угол.
def water_teta(xc, yc, xr, yr):
	y = yc - yr
	x = xc - xr
	if(x==0):
		return(False)
	if(y == 0):
		teta = 0
		return(0)
	teta = y/x
	if( (x > 0) and ( y > 0 )):
		rad = math.atan(teta)
		return(math.degrees(rad))
	if( (x < 0) and (y > 0) ):
		rad = math.atan(teta)
		return(math.degrees(rad) + 180)
	if( (x < 0) and (y < 0) ):
		rad = math.atan(teta)
		return(math.degrees(rad) - 180)
	if( (x > 0) and (y < 0) ):
		rad = math.atan(teta)
		return(math.degrees(rad))



    
# Процедура для 3д вычисляет нужные относительные углы ЛА в пространстве по полученным данным.
# Получает:				1. Угол <alpha_xr> - угол вертикального рассогласования.
#					2. Угол <alpha_xz> - угол горизонтального рассогласования с целью.
#					3. Радиус-вектор. По умолчанию r = 50 м. 
#					4. Предельное отклонение рулей по <xr>.
#					5. Предельное отклонение рулей по <xz>.
#					6. Предельное отклонение рулей по <gamma>.
# Делает:				Вычисляет требуемые углы.
# Возвращает:				Требуемые относительные углы поворота в виде словаря:
#						- <psi> - угол рысканья (горизонтальная плоскость).
#						- <nu>  - угол тангажа (вертикальная плоскость).
#						- <gamma> - угол крена ЛА (вращение вокруг оси движения). 
# Примечание:	система получает два угла и радиус-вектор в сферической связанной СК. Она преобразует эти координаты, выдаёт два угла точно полученные - <psi>, <nu>.
# Они получаются умозрительными заключениями. Также рекомендуется повернуть аппарат на некий угол крена, чтобы увеличть скорость вдоль оси <Z>.
# ДОПОЛНИТЕЛЬНОЕ ПРИМЕЧАНИЕ - процедура анализирует требуемый угол поворота из доступного!
# ДОПОЛНИТЕЛЬНОЕ ПРИМЕЧАНИЕ 2 - процедура вычисляет углы ДИСКРЕТНОГО ПОВОРОТА, т.е. тот угол, на который она может повернуть средствами управления.
def dalpha_3d(alpha_xr, alpha_xz, limxr, limxz, limgamma, r=50):
	
	# Для формовки возврата.
	def get_dret(psi, nu, gamma):
		d = {}
		d['dpsi'] = psi
		d['dnu'] = nu
		d['dgamma'] = gamma
		return(d)
	
	# 1. Если ЛА имеет рассогласование небольшое, которое допустимо поправить лёгким доворотом рулей <psi>, <nu>.
	if((limxr >= abs(alpha_xr))and(limxz >= abs(alpha_xz))):
		#print(limxz)
		dret = get_dret(alpha_xz, alpha_xr, 0)
		return(dret)
		if((alpha_xr >= 0)and(alpha_xz >= 0)):		# Зачем учитывать здесь повторно знак?
			dret = get_dret(alpha_xz, alpha_xr, 0)	# Здесь знак учитывать не надо, т.к. его учли при вычислении рассогласования.
		if((alhpa_xr >= 0)and(alpha_xz < 0)):
			print(-alpha_xz)
			dret = get_dret(alpha_xz, alpha_xr, 0)
		if((alpha_xr < 0)and(alpha_xz < 0)):
			dret = get_dret(alpha_xz, alpha_xr, 0)
			print(-alpha_xz)
		if((alpha_xr < 0)and(alpha_xz >= 0)):
			dret = get_dret(alpha_xz, alpha_xr, 0)
		#print("1. dret = ", dret)
		return(dret)
	# 2. Не хватает по вертикали. 
	elif((limxr < abs(alpha_xr))and(limxz >= abs(alpha_xz))):	# Здесь нужно учесть знак по той причине, что здесь есть поворот на предельный угол.
		if((alpha_xr >= 0)and(alpha_xz >= 0)):
			dret = get_dret(alpha_xz, limxr, 0)
		if((alpha_xr >= 0)and(alpha_xz < 0)):
			dret = get_dret(alpha_xz,limxr, 0)
		if((alpha_xr < 0)and(alpha_xz < 0)):
			dret = get_dret(alpha_xz, -limxr, 0)
		if((alpha_xr < 0)and(alpha_xz >= 0)):
			dret = get_dret(alpha_xz, -limxr, 0)
		return(dret)
	# 3. Не хватает везде. Будет крен. (На самом деле нет, и это здесь показано.)
	elif((limxr < abs(alpha_xr))and(limxz < abs(alpha_xz))):	# И здесь тоже учитывается знак.
		if((alpha_xr >= 0)and(alpha_xz >= 0)):
			dret = get_dret(limxz, limxr, 0)
		if((alpha_xr >= 0)and(alpha_xz < 0)):
			dret = get_dret(-limxz, limxr, 0)
		if((alpha_xr < 0)and(alpha_xz < 0)):
			dret = get_dret(-limxz, -limxr, 0)
		if((alpha_xr < 0)and(limxz >= 0)):
			dret = get_dret(limxz, -limxr, 0)
		return(dret)
	# 4. не хватает горизонта. Крен.
	elif((limxr >= abs(alpha_xr))and(limxz < abs(alpha_xz))):	# Здесь происходит крен.
		if((alpha_xr >= 0)and(alpha_xz >= 0)):
			dret = get_dret(limxz, alpha_xr, 0)
		if((alpha_xr >= 0)and(alpha_xz < 0)):
			dret = get_dret(-limxz, alpha_xr, 0)
		if((alpha_xr < 0)and(alpha_xz < 0)):
			dret = get_dret(-limxz, -alpha_xr, 0)
		if((alpha_xr < 0)and(limxz >= 0)):
			dret = get_dret(limxz, -alpha_xr, 0)
		return(dret)





# Процедура осуществляет теленаведение в одну точку.
# Получает:			1. Список данных о ЛА.
#				2. Список данных о цели.
# Делает:			Рассчитывает требуемые углы полёта ЛА для попадания в заданную точку.
# Возращает:			Словарь из 6-ти углов:
#					1. 'psi' - угол рысканья.
#					2. 'nu'  - угол тангажа.
#					3. 'gamma' - угол крена.
#					4. 'dpsi' - угол поворота на текущем шаге.
#					5. 'dnu' - угол поворота на текущем шаге.
#					6. 'dgamma' - угол поворота на текущем шаге.
# ПРИМЕЧАНИЕ:	углы абсолютные, т.е. как суммарные преподносятся именно углы с учтённым предельным углом поворота в этом ходу.
def f_3d_to_point(lr, lc):
	# 1. Найдём угол абсолютного полёта в системе 360 градусов для каждой из осей.
	xc = lc["x"]
	yc = lc['y']
	zc = lc['z']
	xr = lr['x']
	yr = lr['y']
	zr = lr['z']
	psi = water_teta(xc, zc, xr, zr)	# Этот угол, который ракете нужно держать сейчас, чтобы курс вывел её к цели.
	nu = water_teta(xc, yc, xr, yr)		# Тот же угол.
	
	# 2. Установим точное рассогласование углов.
	delta_nu = nu
	delta_psi = psi
	ret = {}
	ret['psi'] = psi
	ret['nu'] = nu
	ret['gamma'] = 0
	ret['dpsi'] = delta_psi
	ret['dnu'] = delta_nu
	ret['dgamma'] = 0
	return(ret)






# Процедура вычисления углов наведения по 3-ём точкам.
# Получает:		1. Список данных о ЛА.
#			2. Список данных о цели.
# Делает:		Рассчитывает требуемые углы полёта ЛА для попадания по заданной цели.
# Возвращет:		Словарь из 6-ти углов, подробнее см. <f_3d_to_point>
def f_3d_3_point(lr, lc):
	# 1. Забираем координаты ракеты, координаты цели.
	xc = lc['x']
	yc = lc['y']
	zc = lc['z']
	xr = lr['x']
	yr = lr['y']
	zr = lr['z']
	psi = 5




# Процедура вычисления углов наведения по 2-ум тчк, также известен как "собака".
# Получает:		1. Список данных о ЛА.
#			2. Список данных о цели.
# Делает:		Рассчитывает требуемые углы полёта ЛА для попадания по заданной цели.
# Возвращет:		Словарь из 6-ти углов, подробнее см. <f_3d_to_point>
def f_3d_dog(lr, lc):
	# 1. Забираем координаты ракеты, координаты цели.
	xc = lc['x']
	yc = lc['y']
	zc = lc['z']
	xr = lr['x']
	yr = lr['y']
	zr = lr['z']
	
	# 2. Угол требуемого полёта.
	psi = water_teta(xc, zc, xr, zr)
	nu = water_teta(xc, yc, xr, yr)
	
	# 3. Рассогласование углов.
	if(psi != lr['psi']):
		if(psi>0):
			delta_psi = psi - lr['psi']
		else:
			delta_psi = lr['psi'] - psi
	else:
		delta_psi = 0
	if(nu != lr['nu']):
		if(nu>0):
			delta_nu = nu - lr['nu']
		else:
			delta_nu = lr['nu'] - nu
	else:
		delta_nu = 0
	
	# 4. Определить углы полёта за следующий промежуток времени.
	get_alpha = dalpha_3d(delta_nu, delta_psi, lr['limnu'], lr['limpsi'], lr['limgamma'])
	
	# 5. Упаковка ответа.
	ret = {}
	ret['psi'] = lr['psi'] + get_alpha['dpsi']
	ret['nu'] = lr['nu'] + get_alpha['dnu']
	ret['gamma'] = 0
	ret['dpsi'] = get_alpha['dpsi']
	ret['dnu'] = get_alpha['dnu']
	ret['dgamma'] = get_alpha['dgamma']
	return(ret)
	


# Процедура вычисления углов наведения по 2-ум тчк, также известен как "собака".
# Получает:		1. Список данных о ЛА.
#			2. Список данных о цели.
# Делает:		Рассчитывает требуемые углы полёта ЛА для попадания по заданной цели.
# Возвращет:		Словарь из 6-ти углов, подробнее см. <f_3d_to_point>
def f_3d_wolf(lr, lc):
	# 1. Забираем координаты ракеты, координаты цели.
	xc = lc['x']
	yc = lc['y']
	zc = lc['z']
	xr = lr['x']
	yr = lr['y']
	zr = lr['z']
	
	# 1-2. Чтобы взять упреждение, потребуется стрелять с упреждением.
	delta_xc = 1000		# Упреждение по координате.
	# Тупое прибавление к координатам. Логично, если цель летит от ПУ.
	xc += delta_xc
	yc += delta_xc
	zc += delta_xc
	
	
	
	# 2. Угол требуемого полёта.
	psi = water_teta(xc, zc, xr, zr)
	nu = water_teta(xc, yc, xr, yr)
	
	# 3. Рассогласование углов.
	if(psi != lr['psi']):
		if(psi>0):
			delta_psi = psi - lr['psi']
		else:
			delta_psi = lr['psi'] - psi
	else:
		delta_psi = 0
	if(nu != lr['nu']):
		if(nu>0):
			delta_nu = nu - lr['nu']
		else:
			delta_nu = lr['nu'] - nu
	else:
		delta_nu = 0
	
	# 4. Определить углы полёта за следующий промежуток времени.
	get_alpha = dalpha_3d(delta_nu, delta_psi, lr['limnu'], lr['limpsi'], lr['limgamma'])
	
	# 5. Упаковка ответа.
	ret = {}
	ret['psi'] = lr['psi'] + get_alpha['dpsi']
	ret['nu'] = lr['nu'] + get_alpha['dnu']
	ret['gamma'] = 0
	ret['dpsi'] = get_alpha['dpsi']
	ret['dnu'] = get_alpha['dnu']
	ret['dgamma'] = get_alpha['dgamma']
	return(ret)




# Процедура вычисляет абсолютные углы полёта ЛА и является головной процедурой наведения ЛА.
# Получает:			1. Список данных по ЛА:
#					1. x	-	координата ЛА х
#					2. y	-	y
#					3. z	-	z
#					4. Vx	-	скорость ЛА по х
#					5. Vy	-	y
#					6. Vz	-	z
#					7. psi		угол рысканья ЛА
#					8. nu	-	угол тангажа ЛА
#					9. gamma	угол крена ЛА
#					10. limpsi	предельный угол поворота по оси рысканья ЛА за квант времени.
#					11. limnu	предельный угол поворота по оси тангажа ЛА за dt.
#					12. limgamma	предельный угол поворота по оси крена ЛА за dt.
#				2. Список данных по цели:
#					1. x	-	координата цели х
#					2. y	-	у
#					3. z	-	z
#					4. Vx	-	скорость цели по х
#					5. Vy	-	y
#					6. Vz	-	z
#					7. psi		угол рысканья цели
#					8. nu	-	угол тангажа цели
#				3. Флаг решения <letter>. 	Флаг нужен для анализа и передачи в нужный блок.
# Делает:			Анализирует флаг переданного решения и передаёт данные в нужный блок.
# Возвращает:			Углы рысканья, тангажа, крена ЛА на следующий dt в виде словаря:
#					1. 'psi' - угол рысканья.
#					2. 'nu'  - угол тангажа.
#					3. 'gamma' - угол крена.
# Примечание:	Доступные значения флага решения:
#	"to point"	-	наведение в точку.
#	"3 point"	-	3 точки
#	"dog"		-	собачья погоня
#	"wolf"		-	волчья погоня
#	"wisdom"	-	с упреждением.
#	"consciousness"	-	перейти в блок анализа и доверить ему всю черновую работу.
def let_point_path_3d_head(lr, lc, letter):
	if(letter == "to point"):
		# Здесь наведение в точку.
		get = f_3d_to_point(lr, lc)
		return(get)
	if(letter == "3 point"):
		# Здесь наведение по 3-ём точкам.
		get = f_3d_3_point(lr, lc)
		return(get)
	if(letter == "dog"):
		# Здесь наведение собакой.
		get = f_3d_dog(lr, lc)
		return(get)
	if(letter == "wolf"):
		# Здесь наведение волком.
		get = f_3d_wolf(lr, lc)
		return(get)
	if(letter == "wisdom"):
		# Здесь наведение упреждением.
		get = f_3d_wisdom(lr, lc)
		return(get)
	if(letter == "consciousness"):
		# Здесь переход на блок прогноза. 
		# Заглушка для НИРС-2.
		return(0)












# Процедура вычисляет балансировочный угол полёта Р, вычисляет потребный угол поворота по методам наведения, возвращает итоговые углы в этом ходу.
# Получает:		1. P	- тяга в плоскости xy,
#			2. X	- лобовое сопротивление xy,
#			3. Y	- подъёмная сила xy,
#			4. mg	- силу тяжести на ЛА,
#			5. lr	- передаётся в <let_point_thy_path_3d_head>, там и смотри.
#			6. lc	- ...
#			7. letter - ...
#			8. teta_recommend=2	-	Рекомендательный угол атаки.
#				Пояснение: ракета может начать вести активный манёвр, поэтому вычисление угла может быть излишним. Ракета, например, должна будет отработать манёвр, а подбор
#				угла приведёт к получению углов выше 10 градусов. Подбор по методике положительного вертикального ускорения подходит только в том случае, если ракета летит 
#				по +- постоянному углу. В таком случае ракета летит по рекомендательному углу полёта.
#			9. flag_recommend=False	-	Флаг учёта рекомендательного угла.
# Делает:	1. Отправляет данные на расчёт в <let_point_thy_path_3d_head>, оттуда получает углы отклонения ОУ от их нулевого положения.
#		2. Подбирает угол балансировочный полёта в этом ходу.
#		3. Складывает углы и возвращает их.
# Возвращает:		1. teta_xy
#			2. teta_sumxy
#			3. teta_sumxz
#			4. teta_xz
#			5. gamma
#			т.е. конечные углы, которые на слуд.ходу будут закинуты в модель на расчёт.
#			ИЛИ ВОЗВРАЩАЕТ <False>, если был зафиксирован промах.
# ПРИМЕЧАНИЕ: 		смотри описание головной функции наведения.
# ПРИМЕЧАНИЕ_2: 	значения <X>, <Y> в процессе пдбора не меняются, но в реальности меняться будут. Следи за соблюдением инженерной погрешности в 5%!
# ПРИМЕЧАНИЕ_3:		Процедура возвращает <False> если угол атаки непомерно велик (более 10 град).
def guilding_star(P, X, Y, mg, lr, lc, letter, teta_recommend=2, flag_recommend=False):
	
	
	# Функция циклического подбора балансировочного угла.
	def balansir_podbor(P_xy, X_xy, Y_xy, mg, teta_sumxy):
		from trigonometria import sin, cos
		a = -1e10
		teta_xy_sk = -3		# Угол скольжения ракеты.
		while(a <= 0):	# Удельное ускорение всегда больше нуля.
			teta_xy_sk += 0.1	# За итерацию увеличиваем угол.
			teta_sum = teta_sumxy + teta_xy_sk
			# Если угол полёта более 10 град, то ракета не полетит.
			if(teta_xy_sk > 10):
				return False
			teta_sum += teta_xy_sk	# Увеличиваем teta_sumxy.
			a = Y_xy*cos(teta_sumxy - teta_xy_sk) - mg
			#print("balansir_podbor:")
			#print("Y = ", Y)
			#print("\t\ta = %s, teta_sum = %s" %(a, teta_sum))
		return(teta_xy_sk)
	
	
	# 0. Стабилизация углов полёта ракеты.
	# Здесь требуется ввести функцию выборки углов плохого полёта ракеты.
	lr['gamma'] = 0
	
	
	# 1. Вычисление угла наведения.
	# В итоге мы получим список углов именно что требуемого полёта при ПРЕДПОЛОЖЕНИИ полного парирования силы тяжести. 
	# Т.е. углы атаки приняты равными нулю.
	get = let_point_path_3d_head(lr, lc, letter)
	
	ret = {}
	ret['teta_xy'] = 0
	ret['teta_xz'] = 0
	ret['teta_sumxy']  = get['nu']
	ret['teta_sumxz'] = get['psi']
	ret['gamma'] = 0
	
	#print("\tGuilding Star")
	#print("ret = ", ret)
	return(ret)




if __name__ == "__main__":
	P = 100000
	X = 30000
	Y = 30000
	mg = 30000
	lr = {}
	lr['x'] = 0
	lr['y'] = 0
	lr['z'] = 0
	lr['V_x'] = 0
	lr['V_y'] = 0
	lr['V_z'] = 0
	lr['psi'] = 84
	lr['nu'] = 58
	lr['gamma'] = 0
	lr['limpsi'] = 10
	lr['limnu'] = 10
	lr['limgamma'] = 10
	lc = {}
	lc['x'] = 2e4
	lc['y'] = 25e3
	lc['z'] = 150e3
	lc['V_x'] = 0
	lc['V_y'] = 0
	lc['V_z'] = 0
	lc['psi'] = 0
	lc['nu'] = 0
	letter = "to point"
	print(guilding_star(P, X, Y, mg, lr, lc, letter, teta_recommend=2, flag_recommend=True))
