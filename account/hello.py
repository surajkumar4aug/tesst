import datetime
import logging

import constants
import functools
import math
from . import food_recom
from . import cvd_risk
from . import tests_recommend
from . import diet_rules
from . import calc_microbiome
from . import calc_genetics

logger = logging.getLogger(__name__)


def round_helper(num, dec=0):
    num = str(num)[:str(num).find('.')+dec+2]
    if num[-1] >= '5':
        return float(num[:-2-(not dec)]+str(int(num[-2-(not dec)])+1))
    if not len(num[:-1]):
        return 0
    return float(num[:-1])


def cvd_risk_helper(gender, is_smoker, is_diab, age, sbp, cholesterol, ethnicity):
    return cvd_risk.cvd[ethnicity][is_diab][gender][is_smoker]["age"][age]["sbp"][sbp]["cholesterol"][cholesterol]


def secondary_goal_priority_helper(priority_list, count_list):
    high_risk_list = {}
    if len(priority_list) == 5:
        high_risk_list[priority_list[0]] = count_list[0]
        high_risk_list[priority_list[1]] = count_list[1]
        high_risk_list[priority_list[2]] = count_list[2]
        high_risk_list[priority_list[3]] = count_list[3]
        high_risk_list[priority_list[4]] = count_list[4]
    else:
        high_risk_list[priority_list[0]] = count_list[0]
        high_risk_list[priority_list[1]] = count_list[1]
        high_risk_list[priority_list[2]] = count_list[2]
        high_risk_list[priority_list[3]] = count_list[3]
        high_risk_list[priority_list[4]] = count_list[4]
        high_risk_list[priority_list[5]] = count_list[5]
        high_risk_list[priority_list[6]] = count_list[6]
    goal = max(high_risk_list, key=lambda x: high_risk_list[x])
    return goal

def secondary_goal_microbiome_helper(disease_list,data,risk_type):
    risk_list = []
    for i in disease_list:
            for j in range(len(data)):
                for k in range(len(risk_type)):
                    if data[j]["Phenotype"] == i and data[j]["Result"] == risk_type[k]:
                        if risk_type[k] == "Poor":
                            risk_list.append("high_risk")
                        if risk_type[k] == "moderate":
                            risk_list.append("moderate_risk")
                        if risk_type[k] == "Good":
                            risk_list.append("low_risk")
    return risk_list

# def get_food_recommendations_helper(title, gender='', goal=''):
#     food_types = ['veggies', 'carbs', 'fats', 'proteins', 'fruits']
#     result_list = []
#     for food_type in food_types:
#         if goal == "QUE-GST-NED-GH":
#             result_list.append(getattr(food_recom, food_type)[title][gender]['eat'])
#             result_list.append(getattr(food_recom, food_type)[title][gender]['avoid'])
#         else:
#             result_list.append(getattr(food_recom, food_type)[title]['eat'])
#             result_list.append(getattr(food_recom, food_type)[title]['avoid'])
#     return result_list


class Calculate:

    @staticmethod
    def calc_bmi(height, weight):
        bmi = 0
        try:
            bmi = weight/(height/100)**2
        except ZeroDivisionError as ex:
            logger.error(ex)
        return bmi

    @staticmethod
    def bmi_inference(bmi):
        inference = ""
        if bmi <= 18.5:
            inference = "You are under weight*"
        elif 18.5 < bmi <= 24.9:
            inference = "You are of Normal BMI*"

        elif 25.0 <= bmi <= 29.9:
            inference = "You are over weight*"
        elif 30.0 <= bmi <= 39.0:
            inference = "You are Obese*"
        elif bmi >= 40.0:
            inference = "You are severly obese*"
        return inference

    @staticmethod
    def bmi_inference_online(bmi):
        inference = ""
        if bmi <= 18.5:
            inference = 'You are <span style="color: #FDEB11;">Under weight*</span>'
        elif 18.5 < bmi <= 24.9:
            inference = 'Great! You are of <span style="color: #70D44B;">Normal</span> BMI*'

        elif 25.0 <= bmi <= 29.9:
            inference = 'You are <span style="color: #F9B334;">Over weight*</span>'
        elif 30.0 <= bmi <= 39.0:
            inference = 'You are <span style="color: #E94F1B;">Obese*</span>'
        elif bmi >= 40.0:
            inference = 'You are <span style="color: #BF1622;">Severly obese*</span>'
        return inference

    @staticmethod
    def calc_fap(gender, age, bmi):
        if gender == 'F':
            fap = (-44.988) + (0.503 * age) + (10.689 * 1) + (3.172 * bmi) - (0.026 * (bmi)**2) + \
                (0.181 * bmi * 1) - (0.02 * bmi * age) - \
                (0.005 * (bmi)**2 * 1) + (0.00021 * (bmi)**2 * age)
        elif gender == 'M':
            fap = (-44.988) + (0.503 * age) + (10.689 * 0) + (3.172 * bmi) - (0.026 * (bmi)**2) + \
                (0.181 * bmi * 0) - (0.02 * bmi * age) - \
                (0.005 * (bmi)**2 * 0) + (0.00021 * (bmi)**2 * age)

        return fap

    @staticmethod
    def calc_bmr(gender, age, height, weight):
        if gender == 'F':
            bmr = ((10 * weight) + (6.25 * height) - (5 * age)) - 161
        else:
            bmr = ((10 * weight) + (6.25 * height) - (5 * age)) + 5
        return bmr

    @staticmethod
    def calc_rmr(gender, age, height, weight):
        if gender == 'F':
            rmr = 655 + (9.6 * weight) + (1.85 * height) - (4.7 * age)
        else:
            rmr = 66 + (13.7 * weight) + (5 * height) - (6.8 * age)
        return rmr

    @staticmethod
    def calc_spe(bed_time, awake_time, sleep_time, awake_duration, no_of_wakeup):
        if awake_duration == None:
            awake_duration = 0
        if no_of_wakeup == None:
            no_of_wakeup = 0
        if sleep_time == None:
            sleep_time = 0
        else:
            sleep_time = sum(
                i*j for i, j in zip(map(int, sleep_time.split(':')), (60, 1)))
        bed_time = datetime.datetime.strptime(bed_time, '%H:%M')
        awake_time = datetime.datetime.strptime(awake_time, '%H:%M')
        time_diff_in_mins = ((awake_time - bed_time).seconds)//60
        if time_diff_in_mins != 0:
            spe = (time_diff_in_mins - sleep_time -
                   (awake_duration*no_of_wakeup))/time_diff_in_mins
        else:
            spe = 0
        return round(spe*100, 2)

    @staticmethod
    # sleep duraion in hour
    def calc_sleep_duration(bed_time, awake_time, awake_duration, no_of_wakeup, sleep_time):
        if awake_duration == None:
            awake_duration = 0
        if no_of_wakeup == None:
            no_of_wakeup = 0
        if sleep_time == None:
            sleep_time = 0
        else:
            sleep_time = sum(
                i*j for i, j in zip(map(int, sleep_time.split(':')), (60, 1)))
        bed_time = datetime.datetime.strptime(bed_time, '%H:%M')
        awake_time = datetime.datetime.strptime(awake_time, '%H:%M')
        time_diff_in_mins = ((awake_time - bed_time).seconds)//60
        slp = (time_diff_in_mins - sleep_time - (awake_duration*no_of_wakeup))
        return slp/60

    @staticmethod
    def calc_rci(internal_markers):

        return None

    @staticmethod
    def cal_cng_wght(weight, weight_2_months_ago):
        if weight_2_months_ago == None:
            weight_2_months_ago = weight
        return weight - weight_2_months_ago

    @staticmethod
    def cal_snack_wk(frequency_of_snack):
        frequency_of_snack = constants.average_values.get(frequency_of_snack)
        return frequency_of_snack * 7

    @staticmethod
    def cal_snack_sum(frequency_of_fruit, frequency_of_veg, frequency_of_nut, frequency_of_fried_snacks, frequency_of_choco_sweet):
        frequency_of_fruit = constants.average_values.get(frequency_of_fruit)
        frequency_of_veg = constants.average_values.get(frequency_of_veg)
        frequency_of_nut = constants.average_values.get(frequency_of_nut)
        frequency_of_fried_snacks = constants.average_values.get(
            frequency_of_fried_snacks)
        frequency_of_choco_sweet = constants.average_values.get(
            frequency_of_choco_sweet)
        return frequency_of_fruit + frequency_of_veg + frequency_of_nut + frequency_of_fried_snacks + frequency_of_choco_sweet

    @staticmethod
    def cal_snack_qu(frequency_specific, frequency_overall, sum_of_req):
        try:
            return(constants.average_values.get(frequency_specific)/sum_of_req) * (constants.average_values.get(frequency_overall))
            # return (constants.average_values.get(frequency_specific) * constants.average_values.get(frequency_overall))/sum_of_req
        except ZeroDivisionError as ex:
            return 0

    @staticmethod
    def cal_snack_qu_ref(quantity_specific, frequency_specific):
        return constants.average_values.get(quantity_specific) * frequency_specific

    @staticmethod
    def calc_breakfast_wt(breakfast_data):
        breakfast_wt_split = dict()
        breakfast_food_plate_wt = dict()
        breakfast_food_plate_wt_avg = dict()
        weekly_breakfast_avg = dict()
        breakfast_data = dict(
            filter(lambda item: item[1] is not None, breakfast_data.items()))
        # All items internal markers list
        all_items = ["QUE-HEH-MEA-BF-FG-SC", "QUE-HEH-MEA-BF-FG-SB", "QUE-HEH-MEA-BF-FG-EG", "QUE-HEH-MEA-BF-FG-ME", "QUE-HEH-MEA-BF-FG-FS",
                     "QUE-HEH-MEA-BF-FG-NF", "QUE-HEH-MEA-BF-FG-SR", "QUE-HEH-MEA-BF-FG-SL", "QUE-HEH-MEA-BF-FG-SO", "QUE-HEH-MEA-BF-FG-FC", "QUE-HEH-MEA-BF-FG-DA"]
        protein_items = ["QUE-HEH-MEA-BF-FG-SB", "QUE-HEH-MEA-BF-FG-EG", "QUE-HEH-MEA-BF-FG-ME", "QUE-HEH-MEA-BF-FG-FS",
                         "QUE-HEH-MEA-BF-FG-NF"]
        veges_items = ["QUE-HEH-MEA-BF-FG-SR",
                       "QUE-HEH-MEA-BF-FG-SL", "QUE-HEH-MEA-BF-FG-SO"]

        breakfast_satiety = breakfast_data.get("QUE-HEH-SAT-BF", 0)
        breakfast_frequency = breakfast_data.get("QUE-HEH-MEA-BF-FQ", 0)
        breakfast_quantity = constants.meal_quantities.get(
            breakfast_data["gender"], {}).get(breakfast_satiety, 0)

        # sum_of_all_freq = breakfast_data.get("QUE-HEH-MEA-BF-FG-SC",0) + breakfast_data.get("QUE-HEH-MEA-BF-FG-SB",0)+
        sum_of_all_freq = sum(
            filter(None, [breakfast_data.get(k, 0) for k in all_items]))
        sum_of_protein = sum(
            filter(None, [breakfast_data.get(k, 0) for k in protein_items]))
        sum_of_veges = sum(
            filter(None, [breakfast_data.get(k, 0) for k in veges_items]))

        # Breakfast weight split
        for item in all_items:
            try:
                breakfast_wt_split[item] = (
                    breakfast_quantity/sum_of_all_freq)*breakfast_data.get(item, 0)
            except ZeroDivisionError:
                breakfast_wt_split[item] = 0
        sum_of_protein_wt = sum([breakfast_wt_split.get(k, 0)
                                for k in protein_items])
        sum_of_veges_wt = sum([breakfast_wt_split.get(k, 0)
                              for k in veges_items])

        # Food plate based weight
        breakfast_food_plate_wt["QUE-HEH-MEA-BF-FP-GR"] = (
            breakfast_quantity/100)*breakfast_data.get("QUE-HEH-MEA-BF-FP-GR", 0)
        breakfast_food_plate_wt["QUE-HEH-MEA-BF-FP-PR"] = (
            breakfast_quantity/100)*breakfast_data.get("QUE-HEH-MEA-BF-FP-PR", 0)
        breakfast_food_plate_wt["QUE-HEH-MEA-BF-FP-VG"] = (
            breakfast_quantity/100)*breakfast_data.get("QUE-HEH-MEA-BF-FP-VG", 0)
        breakfast_food_plate_wt["QUE-HEH-MEA-BF-FP-FR"] = (
            breakfast_quantity/100)*breakfast_data.get("QUE-HEH-MEA-BF-FP-FR", 0)
        breakfast_food_plate_wt["QUE-HEH-MEA-BF-FP-DR"] = (
            breakfast_quantity/100)*breakfast_data.get("QUE-HEH-MEA-BF-FP-DR", 0)

        # Food plate based weight average
        breakfast_food_plate_wt_avg["QUE-HEH-MEA-BF-FG-SC"] = breakfast_food_plate_wt.get(
            "QUE-HEH-MEA-BF-FP-GR", 0)

        for pr_item in protein_items:
            try:
                breakfast_food_plate_wt_avg[pr_item] = (breakfast_food_plate_wt.get(
                    "QUE-HEH-MEA-BF-FP-PR", 0)/sum_of_protein_wt)*breakfast_wt_split.get(pr_item, 0)
            except ZeroDivisionError:
                breakfast_food_plate_wt_avg[pr_item] = 0
        for veg_item in veges_items:
            try:
                breakfast_food_plate_wt_avg[veg_item] = (breakfast_food_plate_wt.get(
                    "QUE-HEH-MEA-BF-FP-VG", 0)/sum_of_veges_wt)*breakfast_wt_split.get(veg_item, 0)
            except ZeroDivisionError:
                breakfast_food_plate_wt_avg[veg_item] = 0

        breakfast_food_plate_wt_avg["QUE-HEH-MEA-BF-FG-FC"] = breakfast_food_plate_wt.get(
            "QUE-HEH-MEA-BF-FP-FR", 0)
        breakfast_food_plate_wt_avg["QUE-HEH-MEA-BF-FG-DA"] = breakfast_food_plate_wt.get(
            "QUE-HEH-MEA-BF-FP-DR", 0)
        # Weekly Food plate based weight average
        for item in all_items:
            weekly_breakfast_avg[item+"-AVG"] = (
                breakfast_food_plate_wt_avg.get(item, 0)*breakfast_frequency)/7
        return weekly_breakfast_avg

    @staticmethod
    def calc_lunch_wt(lunch_data):
        lunch_wt_split = dict()
        lunch_food_plate_wt = dict()
        lunch_food_plate_wt_avg = dict()
        weekly_lunch_avg = dict()

        lunch_data = dict(
            filter(lambda item: item[1] is not None, lunch_data.items()))

        all_items = ["QUE-HEH-MEA-LN-FG-CE", "QUE-HEH-MEA-LN-FG-BE", "QUE-HEH-MEA-LN-FG-EG", "QUE-HEH-MEA-LN-FG-ME", "QUE-HEH-MEA-LN-FG-FI",
                     "QUE-HEH-MEA-LN-FG-NF", "QUE-HEH-MEA-LN-FG-RO", "QUE-HEH-MEA-LN-FG-GR", "QUE-HEH-MEA-LN-FG-OT", "QUE-HEH-MEA-LN-FG-FC", "QUE-HEH-MEA-LN-FG-DA"]
        protein_items = ["QUE-HEH-MEA-LN-FG-BE", "QUE-HEH-MEA-LN-FG-EG", "QUE-HEH-MEA-LN-FG-ME", "QUE-HEH-MEA-LN-FG-FI",
                         "QUE-HEH-MEA-LN-FG-NF"]
        veges_items = ["QUE-HEH-MEA-LN-FG-RO",
                       "QUE-HEH-MEA-LN-FG-GR", "QUE-HEH-MEA-LN-FG-OT"]

        lunch_satiety = lunch_data.get("QUE-HEH-SAT-LN", 0)
        lunch_frequency = lunch_data.get("QUE-HEH-MEA-LN-FQ", 0)
        lunch_quantity = constants.meal_quantities.get(
            lunch_data["gender"], {}).get(lunch_satiety, 0)

        # sum_of_all_freq = lunch_data.get("QUE-HEH-MEA-BF-FG-SC",0) + lunch_data.get("QUE-HEH-MEA-LN-FG-BE",0)+
        sum_of_all_freq = sum(
            filter(None, [lunch_data.get(k, 0) for k in all_items]))
        sum_of_protein_freq = sum(
            filter(None, [lunch_data.get(k, 0) for k in protein_items]))
        sum_of_veges_freq = sum(
            filter(None, [lunch_data.get(k, 0) for k in veges_items]))
        # lunch weight split
        for item in all_items:
            try:
                lunch_wt_split[item] = (
                    lunch_quantity/sum_of_all_freq)*lunch_data.get(item, 0)
            except ZeroDivisionError:
                lunch_wt_split[item] = 0
        sum_of_protein_wt = sum([lunch_wt_split.get(k, 0)
                                for k in protein_items])
        sum_of_veges_wt = sum([lunch_wt_split.get(k, 0) for k in veges_items])

        # Food plate based weight

        lunch_food_plate_wt["QUE-HEH-MEA-LN-FP-GR"] = (
            lunch_quantity/100)*lunch_data.get("QUE-HEH-MEA-LN-FP-GR", 0)
        lunch_food_plate_wt["QUE-HEH-MEA-LN-FP-PR"] = (
            lunch_quantity/100)*lunch_data.get("QUE-HEH-MEA-LN-FP-PR", 0)
        lunch_food_plate_wt["QUE-HEH-MEA-LN-FP-VG"] = (
            lunch_quantity/100)*lunch_data.get("QUE-HEH-MEA-LN-FP-VG", 0)
        lunch_food_plate_wt["QUE-HEH-MEA-LN-FP-FR"] = (
            lunch_quantity/100)*lunch_data.get("QUE-HEH-MEA-LN-FP-FR", 0)
        lunch_food_plate_wt["QUE-HEH-MEA-LN-FP-DR"] = (
            lunch_quantity/100)*lunch_data.get("QUE-HEH-MEA-LN-FP-DR", 0)
        # Food plate based weight average
        lunch_food_plate_wt_avg["QUE-HEH-MEA-LN-FG-CE"] = lunch_food_plate_wt.get(
            "QUE-HEH-MEA-LN-FP-GR", 0)

        for pr_item in protein_items:
            try:
                lunch_food_plate_wt_avg[pr_item] = (lunch_food_plate_wt.get(
                    "QUE-HEH-MEA-LN-FP-PR", 0)/sum_of_protein_wt)*lunch_wt_split.get(pr_item, 0)
            except ZeroDivisionError:
                lunch_food_plate_wt_avg[pr_item] = 0
        for veg_item in veges_items:
            try:
                lunch_food_plate_wt_avg[veg_item] = (lunch_food_plate_wt.get(
                    "QUE-HEH-MEA-LN-FP-VG", 0)/sum_of_veges_wt)*lunch_wt_split.get(veg_item, 0)
            except ZeroDivisionError:
                lunch_food_plate_wt_avg[veg_item] = 0

        lunch_food_plate_wt_avg["QUE-HEH-MEA-LN-FG-FC"] = lunch_food_plate_wt.get(
            "QUE-HEH-MEA-LN-FP-FR", 0)
        lunch_food_plate_wt_avg["QUE-HEH-MEA-LN-FG-DA"] = lunch_food_plate_wt.get(
            "QUE-HEH-MEA-LN-FP-DR", 0)
        # Weekly Food plate based weight average
        for item in all_items:
            weekly_lunch_avg[item+"-AVG"] = (
                lunch_food_plate_wt_avg.get(item, 0)*lunch_frequency)/7
        return weekly_lunch_avg

    @staticmethod
    def calc_dinner_wt(dinner_data):
        dinner_wt_split = dict()
        dinner_food_plate_wt = dict()
        dinner_food_plate_wt_avg = dict()
        weekly_dinner_avg = dict()
        dinner_data = dict(
            filter(lambda item: item[1] is not None, dinner_data.items()))

        all_items = ["QUE-HEH-MEA-DN-FG-CE", "QUE-HEH-MEA-DN-FG-BE", "QUE-HEH-MEA-DN-FG-EG", "QUE-HEH-MEA-DN-FG-ME", "QUE-HEH-MEA-DN-FG-FI",
                     "QUE-HEH-MEA-DN-FG-NF", "QUE-HEH-MEA-DN-FG-RO", "QUE-HEH-MEA-DN-FG-GR", "QUE-HEH-MEA-DN-FG-OT", "QUE-HEH-MEA-DN-FG-FC", "QUE-HEH-MEA-DN-FG-DA"]
        protein_items = ["QUE-HEH-MEA-DN-FG-BE", "QUE-HEH-MEA-DN-FG-EG", "QUE-HEH-MEA-DN-FG-ME", "QUE-HEH-MEA-DN-FG-FI",
                         "QUE-HEH-MEA-DN-FG-NF"]
        veges_items = ["QUE-HEH-MEA-DN-FG-RO",
                       "QUE-HEH-MEA-DN-FG-GR", "QUE-HEH-MEA-DN-FG-OT"]

        dinner_satiety = dinner_data.get("QUE-HEH-SAT-DN", 0)
        dinner_frequency = dinner_data.get("QUE-HEH-MEA-DN-FQ", 0)
        dinner_quantity = constants.meal_quantities.get(
            dinner_data["gender"], {}).get(dinner_satiety, 0)

        # sum_of_all_freq = dinner_data.get("QUE-HEH-MEA-BF-FG-SC",0) + dinner_data.get("QUE-HEH-MEA-DN-FG-BE",0)+
        sum_of_all_freq = sum(
            filter(None, [dinner_data.get(k, 0) for k in all_items]))
        sum_of_protein_freq = sum(
            filter(None, [dinner_data.get(k, 0) for k in protein_items]))
        sum_of_veges_freq = sum(
            filter(None, [dinner_data.get(k, 0) for k in veges_items]))

        # dinner weight split
        for item in all_items:
            try:
                dinner_wt_split[item] = (
                    dinner_quantity/sum_of_all_freq)*dinner_data.get(item, 0)
            except ZeroDivisionError:
                dinner_wt_split[item] = 0
        sum_of_protein_wt = sum([dinner_wt_split.get(k, 0)
                                for k in protein_items])
        sum_of_veges_wt = sum([dinner_wt_split.get(k, 0) for k in veges_items])

        # Food plate based weight
        dinner_food_plate_wt["QUE-HEH-MEA-DN-FP-GR"] = (
            dinner_quantity/100)*dinner_data.get("QUE-HEH-MEA-DN-FP-GR", 0)
        dinner_food_plate_wt["QUE-HEH-MEA-DN-FP-PR"] = (
            dinner_quantity/100)*dinner_data.get("QUE-HEH-MEA-DN-FP-PR", 0)
        dinner_food_plate_wt["QUE-HEH-MEA-DN-FP-VG"] = (
            dinner_quantity/100)*dinner_data.get("QUE-HEH-MEA-DN-FP-VG", 0)
        dinner_food_plate_wt["QUE-HEH-MEA-DN-FP-FR"] = (
            dinner_quantity/100)*dinner_data.get("QUE-HEH-MEA-DN-FP-FR", 0)
        dinner_food_plate_wt["QUE-HEH-MEA-DN-FP-DR"] = (
            dinner_quantity/100)*dinner_data.get("QUE-HEH-MEA-DN-FP-DR", 0)

        # Food plate based weight average
        dinner_food_plate_wt_avg["QUE-HEH-MEA-DN-FG-CE"] = dinner_food_plate_wt.get(
            "QUE-HEH-MEA-DN-FP-GR", 0)

        for pr_item in protein_items:
            try:

                dinner_food_plate_wt_avg[pr_item] = (dinner_food_plate_wt.get(
                    "QUE-HEH-MEA-DN-FP-PR", 0)/sum_of_protein_wt)*dinner_wt_split.get(pr_item, 0)
            except ZeroDivisionError:
                dinner_food_plate_wt_avg[pr_item] = 0
        for veg_item in veges_items:
            try:
                dinner_food_plate_wt_avg[veg_item] = (dinner_food_plate_wt.get(
                    "QUE-HEH-MEA-DN-FP-VG", 0)/sum_of_veges_wt)*dinner_wt_split.get(veg_item, 0)
            except ZeroDivisionError:
                dinner_food_plate_wt_avg[veg_item] = 0

        dinner_food_plate_wt_avg["QUE-HEH-MEA-DN-FG-FC"] = dinner_food_plate_wt.get(
            "QUE-HEH-MEA-DN-FP-FR", 0)
        dinner_food_plate_wt_avg["QUE-HEH-MEA-DN-FG-DA"] = dinner_food_plate_wt.get(
            "QUE-HEH-MEA-DN-FP-DR", 0)

        # Weekly Food plate based weight average
        for item in all_items:
            weekly_dinner_avg[item+"-AVG"] = (
                dinner_food_plate_wt_avg.get(item, 0)*dinner_frequency)/7

        return weekly_dinner_avg

    def cal_food_ctgry(quantity_specific, daily_intake, adjustment_factor):
        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Energy"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        if daily_intake is not None:
            energy = round((((dnt_sum)/dnt_len) / 100) *
                           daily_intake if (dnt_sum) != 0 else 0, 1)
        else:
            energy = 0   
        adjusted_energy = energy + \
            (energy*adjustment_factor) if dnt_sum != 0 else 0
        adjusted_daily_intake = adjusted_energy * \
            (daily_intake/energy) if energy != 0 else 0
        sec_energy = ((dnt_sum/dnt_len) / 100) * adjusted_daily_intake if (
            sum(dnt_values)) != 0 else 0 if sum(dnt_values) != 0 else 0
        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "CHO"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        cho = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Ptn"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        ptn = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Fat"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        fat = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Fiber"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        fiber = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0
        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Vit A"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        vita = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0
        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Vit B1"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        vitb1 = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Vit B2"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        vitb2 = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Vit B3"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        vitb3 = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Vit B5"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        vitb5 = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Vit B6"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        vitb6 = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Vit B7-Biotin"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        vitb7 = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Vit B9-Folate"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        vitb9 = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Vit B12"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        vitb12 = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Vit C"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        vitc = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Vit D"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        vitd = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Vit E"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        vite = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0
        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Vit K"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        vitk = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Iron"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        iron = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Zinc"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        zinc = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Magnesium"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        magnesium = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Calcium"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        calcium = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Selenium"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        selenium = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Copper"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        copper = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Manganese"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        manganes = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Chromium"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        chromium = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Phosphorus"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        phosporous = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Potassium"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        potassium = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Sodium"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        sodium = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Saturated fats"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        saturated = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Trans fats"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        trans = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "MUFA"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        mufa = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "PUFA"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        pufa = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Total Sugars"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        totalsugar = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Sugars added"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        sugarad = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Cholesterol"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        cholester = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Water"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        water = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Choline"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        choline = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "caffeine"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        caffeine = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0

        dnt_values = [y for i in quantity_specific for y in constants.dnt.get(i)[
            "Alcohol"]]
        dnt_sum = sum(dnt_values)
        dnt_len = len(dnt_values) if sum != 0 else 1
        alcohol = ((dnt_sum/dnt_len) / 100) * \
            adjusted_daily_intake if dnt_sum != 0 else 0
        return adjusted_energy, adjusted_daily_intake, sec_energy, cho, ptn, fat, fiber, vita, vitb1, vitb2, vitb3, vitb5, vitb6, vitb7, vitb9, vitb12, vitc, vitd, vite, vitk, iron, zinc, magnesium, calcium, selenium, copper, manganes, chromium, phosporous, potassium, sodium, saturated, trans, mufa, pufa, totalsugar, sugarad, cholester, water, choline, caffeine, alcohol

    @staticmethod
    def cal_food_ctgry_sum(specific_values_category):
        return functools.reduce(lambda a, b: a+b, specific_values_category)

    @staticmethod
    def calc_pal(activity_data, sleep_data):
        pal_dict = dict()
        sum_product_of_exc = 0
        total_time_in_day = 1440  # In mins
        exercise = activity_data.get("QUE-HEH-EAL-VO-M", 0)
        exercise_par = activity_data.get("QUE-HEH-EAL-PAR", 0)
        household = activity_data.get("QUE-HEH-HPAL-FQ-H", 0)
        household_par = activity_data.get("QUE-HEH-HPAL-FQ", 0)
        occupation = activity_data.get("QUE-HEH-PAL-FQ-M", 0)
        occupation_par = activity_data.get("QUE-HEH-PAL-PAR", 0)
        sleep = sleep_data.get("QUE-GST-SLP-TS", 0)
        sleep_par = 1
        others_without_par = total_time_in_day - \
            (sum(exercise) + household + occupation + sleep)
        others_with_par = others_without_par * 1.2

        for index in range(len(exercise)):
            sum_product_of_exc += exercise[index] * exercise_par[index]
        sum_product_of_all = sum_product_of_exc + (household * household_par) + (
            occupation * occupation_par) + (sleep * sleep_par) + (others_with_par)

        pal = sum_product_of_all / total_time_in_day
        if sum_product_of_all == 0:
            pal_dict["exercise_exp"] = 0
            pal_dict["sleep_exp"] = 0
            pal_dict["occupation_exp"] = 0
            pal_dict["household_exp"] = 0
            pal_dict["other_exp"] = 0
            # pal_dict["pal"] = round(pal,2)
            pal_dict["pal"] = pal

        else:
            # pal_dict["pal"] = round(pal,2)
            pal_dict["pal"] = pal

            exercise_exp = (sum_product_of_exc/sum_product_of_all)*100
            sleep_exp = (sleep/sum_product_of_all)*100
            occupation_exp = ((occupation * occupation_par) /
                              sum_product_of_all)*100
            household_exp = ((household * household_par) /
                             sum_product_of_all)*100
            other_exp = (others_with_par/sum_product_of_all)*100
            pal_dict["exercise_exp"] = round(exercise_exp, 0)
            pal_dict["sleep_exp"] = round(sleep_exp, 0)
            pal_dict["occupation_exp"] = round(occupation_exp, 0)
            pal_dict["household_exp"] = round(household_exp, 0)
            pal_dict["other_exp"] = round(other_exp, 0)
        return pal_dict

    @staticmethod
    def calc_tdee(weight_change):
        # weight change is required from hpe.
        # weight_change = 2
        wt = weight_change/2
        diff_per_day = (wt*7000)/30
        return diff_per_day

    @staticmethod
    def calc_adj_tdee(bmr, pal, weight_change_reason, tdee):
        if weight_change_reason !=None:
            return round((bmr*pal) + tdee, 0)
        else:
            return round((bmr*pal), 0)

    @staticmethod
    def calc_exc_def_cal(adj_tdee, ei):
        return adj_tdee - ei

    @staticmethod
    def calc_adj_factor(exc_def_cal, ei):
        try:
            return exc_def_cal/ei
        except ZeroDivisionError:
            return 0

    @staticmethod
    def calc_activity_req(weight, exc_def_cal):
        cal_burn_light = (2*3.5*weight)/200
        cal_burn_mod = (4*3.5*weight)/200
        cal_burn_vig = (8*3.5*weight)/200

        try:
            min_of_exc_req_light = exc_def_cal/cal_burn_light
            min_of_exc_req_mod = exc_def_cal/cal_burn_mod
            min_of_exc_req_vig = exc_def_cal/cal_burn_vig
        except ZeroDivisionError:
            min_of_exc_req_light = 0
            min_of_exc_req_mod = 0
            min_of_exc_req_vig = 0
        return {"NPE-RUE-CALB-L": cal_burn_light, "NPE-RUE-CALB-M": cal_burn_mod, "NPE-RUE-CALB-V": cal_burn_vig,
                "NPE-RUE-EXCR-L": min_of_exc_req_light, "NPE-RUE-EXCR-M": min_of_exc_req_mod, "NPE-RUE-EXCR-V": min_of_exc_req_vig}

    @staticmethod
    def cal_food_ctgry_per(sum_value, admr):
        return sum_value*admr

    @staticmethod
    def cal_food_ctgry_sum_per(specific_values_category, all_values):
        try:
            functools.reduce(lambda a, b: a+b, all_values)*100, 2 != 0
            return round((specific_values_category/functools.reduce(lambda a, b: a+b, all_values))*100, 2)
        except:
            return 0

    @staticmethod
    def cal_food_ctgry_src(sum, specific_value):
        return (specific_value/sum)*100 if sum != 0 else 0

    @staticmethod
    def cal_food_per_veg_fru(sum, specific_value):
        return (functools.reduce(lambda a, b: a+b, specific_value)/sum)*100 if sum != 0 else 0

    @staticmethod
    def cal_food_deficit(estimated_intake, specific_category, gender, age, pregnant, lactating, smoker):
        gender = "female" if gender is 'F' else "male"
        rda_constants = constants.rda.get(gender).get(specific_category, {})
        rda_constant = Calculate.get_rda(
            age, rda_constants, pregnant, lactating, smoker,specific_category)
        # deficit = estimated_intake - rda_constant if rda_constant != "" else ""
        try:
            rda_met = (estimated_intake / rda_constant)*100
        except ZeroDivisionError:
            rda_met = 0
        deficit = 100 - rda_met if rda_met < 100 else 0
        deficit = abs(round(deficit, 2)) if isinstance(
            deficit, float) else deficit
        rda_met = round(rda_met, 2) if isinstance(rda_met, float) else rda_met
        status_color = "red" if rda_met < 100 else "green"
        return deficit, rda_met, status_color

    @staticmethod
    def get_rda(age, rda_constants, pregnant, lactating, smoker,specific_category):
        above_51 = ["Vit A","Vit B1","Vit B2","Vit B6","Iron","Selenium","Chromium","Manganese","Potassium"]
        above_19 = ["Vit B3","Vit B5","Vit B7-Biotin","Vit B9-Folate","Vit B12","Vit C","Vit K","Zinc","Copper","Phosphorus","Choline"]
        above_71 = ["Vit D","Calcium"]
        above_14 = ["Vit E"]
        above_18 = ["Sodium","Cholesterol","Water"]
        caffine  = ["caffeine"]
        alcohol  = ["Alcohol"]
        magnesium = ["Magnesium"]
        rda_constant = 0
        if specific_category in magnesium:
            if 0 < age <= 0.6:
                rda_constant = rda_constants.get("0-0.6").get("value")
            elif 0.6 < age < 1:
                rda_constant = rda_constants.get("0.7-1").get("value")
            elif 1 <= age <= 3:
                rda_constant = rda_constants.get("1-3").get("value")
            elif 4 < age <= 8:
                rda_constant = rda_constants.get("4-8").get("value")
            elif 8 < age <= 13:
                rda_constant = rda_constants.get("9-13").get("value")
            elif 14 <= age <= 18:
                rda_constant = rda_constants.get("14-18").get("value")
            elif 19 <= age <= 30:
                rda_constant = rda_constants.get("19-30").get("value")
            elif 31 <= age <= 50:
                rda_constant = rda_constants.get("31-50").get("value")
            elif age >= 51:
                rda_constant = rda_constants.get("51+").get("value") 
        elif specific_category in above_51:
            if 0 < age <= 0.6:
                rda_constant = rda_constants.get("0-0.6").get("value")
            elif 0.6 < age < 1:
                rda_constant = rda_constants.get("0.7-1").get("value")
            elif 1 <= age <= 3:
                rda_constant = rda_constants.get("1-3").get("value")
            elif 4 < age <= 8:
                rda_constant = rda_constants.get("4-8").get("value")
            elif 8 < age <= 13:
                rda_constant = rda_constants.get("9-13").get("value")
            elif 13 < age <= 18:
                rda_constant = rda_constants.get("14-18").get("value")
            elif 19 <= age <= 50:
                rda_constant = rda_constants.get("19-50").get("value")
            elif age >= 51:
                rda_constant = rda_constants.get("51+").get("value")
        elif specific_category in above_19:
            if 0 < age <= 0.6:
                rda_constant = rda_constants.get("0-0.6").get("value")
            elif 0.6 < age < 1:
                rda_constant = rda_constants.get("0.7-1").get("value")
            elif 1 <= age <= 3:
                rda_constant = rda_constants.get("1-3").get("value")
            elif 4 < age <= 8:
                rda_constant = rda_constants.get("4-8").get("value")
            elif 8 < age <= 13:
                rda_constant = rda_constants.get("9-13").get("value")
            elif 14 < age <= 18:
                rda_constant = rda_constants.get("14-18").get("value")
            elif age > 19:
                rda_constant = rda_constants.get("19+").get("value")
        elif specific_category in above_71:
            if 0 < age <= 0.6:
                rda_constant = rda_constants.get("0-0.6").get("value")
            elif 0.6 < age < 1:
                rda_constant = rda_constants.get("0.7-1").get("value")
            elif 1 <= age <= 3:
                rda_constant = rda_constants.get("1-3").get("value")
            elif 4 < age <= 8:
                print(age)
                #rda_constant = rda_constants.get("4-8").get("value")
            elif 8 < age <= 13:
                rda_constant = rda_constants.get("9-13").get("value")
            elif 14 < age <= 18:
                rda_constant = rda_constants.get("14-18").get("value")
            elif 19 < age <= 50:
                rda_constant = rda_constants.get("19-50").get("value")
            elif 51 <= age <= 70:
                rda_constant = rda_constants.get("51-70").get("value")
            elif age > 71:
                rda_constant = rda_constants.get("71+").get("value")
        elif specific_category in above_14:
            if 0 < age <= 0.6:
                rda_constant = rda_constants.get("0-0.6").get("value")
            elif 0.6 < age < 1:
                rda_constant = rda_constants.get("0.7-1").get("value")
            elif 1 <= age <= 3:
                rda_constant = rda_constants.get("1-3").get("value")
            elif 4 < age <= 8:
                rda_constant = rda_constants.get("4-8").get("value")
            elif 8 < age <= 13:
                rda_constant = rda_constants.get("9-13").get("value")
            elif age >14:
                rda_constant = rda_constants.get("14+").get("value")

        elif specific_category in above_18:
            if age >= 18:
                rda_constant = rda_constants.get("18+").get("value")
        elif specific_category in caffine:
            if age >= 18:
                rda_constant = rda_constants.get("18+").get("value")
        elif specific_category in alcohol:
            if age >= 19:   
                rda_constant = rda_constants.get("19+").get("value") 

        if rda_constant != 0:
            # if pregnant:
            #     rda_constant += rda_constants.get("pregnant")
            # if lactating:
            #     rda_constant += rda_constants.get("lactating")
            if smoker:
                rda_constant += 35

        return rda_constant

    @staticmethod
    def cal_total_intake(breakfast_intake, lunch_intake, dinner_intake):
        cereal_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-SC-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-CE-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-CE-AVG", 0)
        beans_and_leg_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-SB-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-BE-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-BE-AVG", 0)
        
        egg_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-EG-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-EG-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-EG-AVG", 0)

        meat_poultry_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-ME-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-ME-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-ME-AVG", 0)

        fish_seafood_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-FS-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-FI-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-FI-AVG", 0)

        nuts_dry_fruits_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-NF-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-NF-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-NF-AVG", 0)

        root_starchy_veg_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-SR-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-RO-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-RO-AVG", 0)

        green_leafy_veg_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-SL-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-GR-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-GR-AVG", 0)

        other_veg_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-SO-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-OT-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-OT-AVG", 0)

        fresh_fruits_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-FC-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-FC-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-FC-AVG", 0)

        dairy_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-DA-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-DA-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-DA-AVG", 0)    

        intake_total = {"QUE-HEH-MEA-CE-TOT": cereal_total, "QUE-HEH-MEA-BE-TOT": beans_and_leg_total, "QUE-HEH-MEA-EG-TOT": egg_total,
                        "QUE-HEH-MEA-ME-TOT": meat_poultry_total, "QUE-HEH-MEA-FI-TOT": fish_seafood_total, "QUE-HEH-MEA-NF-TOT": nuts_dry_fruits_total,
                        "QUE-HEH-MEA-RO-TOT": root_starchy_veg_total, "QUE-HEH-MEA-GR-TOT": green_leafy_veg_total, "QUE-HEH-MEA-OT-TOT": other_veg_total,
                        "QUE-HEH-MEA-FC-TOT": fresh_fruits_total, "QUE-HEH-MEA-DA-TOT": dairy_total}
        return intake_total

    @staticmethod
    def cal_total_veg_intake(breakfast_intake, lunch_intake, dinner_intake):

        root_starchy_veg_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-SR-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-RO-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-RO-AVG", 0)

        green_leafy_veg_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-SL-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-GR-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-GR-AVG", 0)

        other_veg_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-SO-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-OT-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-OT-AVG", 0)

        fresh_fruits_total = breakfast_intake.get("QUE-HEH-MEA-BF-FG-FC-AVG", 0) + lunch_intake.get("QUE-HEH-MEA-LN-FG-FC-AVG", 0) + \
            dinner_intake.get("QUE-HEH-MEA-DN-FG-FC-AVG", 0)

        intake_total = {
            "QUE-HEH-MEA-RO-TOT": root_starchy_veg_total, "QUE-HEH-MEA-GR-TOT": green_leafy_veg_total, "QUE-HEH-MEA-OT-TOT": other_veg_total,
            "QUE-HEH-MEA-FC-TOT": fresh_fruits_total}

        return intake_total

    @staticmethod
    def cal_food_energy(quantity_specific, daily_intake):
        try:    
            energy_dnt_sum = 0
            dnt_len_sum = 0
            for type in quantity_specific:
                energy_dnts = constants.dnt.get(type, "").get("Energy")
                energy_dnt_sum += sum(energy_dnts)
                dnt_len_sum += len(energy_dnts)
            energy = ((energy_dnt_sum/dnt_len_sum)/100)*daily_intake
            # energy =  ( ((sum(constants.dnt.get(i)["Energy"] for i in quantity_specific))/len(quantity_specific)) / 100 ) * daily_intake if (sum(constants.dnt.get(i)["Energy"] for i in quantity_specific)) != 0 else 0
            return energy
        except:
            return 0

    @staticmethod
    def get_snack(snack_list):
        if snack_list.get("QUE-HEH-MEA-SN-FQ-1") == True:
            return "QUE-HEH-MEA-SN-FQ-1"
        elif snack_list.get("QUE-HEH-MEA-SN-FQ-2") == True:
            return "QUE-HEH-MEA-SN-FQ-2"
        elif snack_list.get("QUE-HEH-MEA-SN-FQ-3") == True:
            return "QUE-HEH-MEA-SN-FQ-3"
        elif snack_list.get("QUE-HEH-MEA-SN-FQ-45") == True:
            return "QUE-HEH-MEA-SN-FQ-45"
        else:
            return "QUE-HEH-MEA-SN-FQ-NV"

    @staticmethod
    def get_diabets_risk(gender, meds, age, bmi, smoke, diabetes):
        a = -6.322
        gen_value = 0 if gender == "M" else -0.879
        medication = meds.get("QUE-GST-MED-AM")
        ah_med_value = 1.222 if medication == True else 0
        steroids = meds.get("QUE-GST-MED-S")
        ster_value = 2.191 if steroids == True else 0
        age_value = 0 if age < 25 else age * 0.063

        if bmi < 25:
            bmi_value = 0
        elif 25 <= bmi <= 27.5:
            bmi_value = 0.699
        elif 27.5 <= bmi <= 30:
            bmi_value = 1.97
        elif bmi >= 30:
            bmi_value = 2.518

        if smoke.get('QUE-GST-SMK-Y') == True:
            smoke_value = 0.855
        elif smoke.get('QUE-GST-SMK-X') == True:
            smoke_value = -0.218
        elif smoke.get('QUE-GST-SMK-N') == True:
            smoke_value = 0
        else:
            smoke_value = 0

        if diabetes.get("QUE-GST-DIA-Y-OR") == True:
            dia_value = 0.728
        elif diabetes.get("QUE-GST-DIA-Y-AN") == True:
            dia_value = 0.753
        elif diabetes.get("QUE-GST-DIA-N") == True:
            dia_value = 0
        else:
            dia_value = 0
        exp = math.exp(-(a + (gen_value) + (ah_med_value) + (ster_value) +
                       (age_value) + (bmi_value) + (dia_value) + (smoke_value)))
        diab_risk = 1/(1+exp)
        return round(diab_risk*100, 1)

    @staticmethod
    def lung_disease_risk(smoke_data):
        try:
            smoke_per_day = int(smoke_data.get("QUE-GST-SMK-STK-S-FQ"))
            smoke_years = smoke_data.get("QUE-GST-SMK-STK-Y")
            pack_years = smoke_per_day / 20 * smoke_years
            pack_years = round_helper(pack_years)
        except:
            pack_years = 0

        return round(pack_years)

    # @staticmethod
    # def get_food_recommendations(consumer_internal_markers, gender):
    #     goal = [k for k, v in consumer_internal_markers.items() if v == True]
    #     selected_gender = 'female' if gender == 'F' else 'male'
    #     goal_details_dict = {
    #         "QUE-GST-NED-GH": get_food_recommendations_helper('general_health', selected_gender, goal[0]),
    #         "QUE-GST-NED-BR": get_food_recommendations_helper('brain_health'),
    #         "QUE-GST-NED-GF": get_food_recommendations_helper('gastro_fitness'),
    #         "QUE-GST-NED-PCN": get_food_recommendations_helper('post_covid'),
    #         "QUE-GST-NED-PE": get_food_recommendations_helper('performance_and_endurance'),
    #         "QUE-GST-NED-SK": get_food_recommendations_helper('skin'),
    #         "QUE-GST-NED-SP": get_food_recommendations_helper('sleep'),
    #         "QUE-GST-NED-SW": get_food_recommendations_helper('shiftworker'),
    #         "QUE-GST-NED-WL": get_food_recommendations_helper('weight_reduction')
    #     }
    #     return goal_details_dict[goal[0]]
        # return veggies_to_eat, veggies_to_avoid, proteins_to_eat, proteins_to_avoid, fats_to_eat, fats_to_avoid, carbs_to_eat, carbs_to_avoid, fruits_to_eat, fruits_to_avoid
    @staticmethod
    def get_food_recommendations(consumer_internal_markers, gender):
        goal = [k for k, v in consumer_internal_markers.items() if v == True]
        if goal == ["QUE-GST-NED-GH"]:
            if gender == "F":
                veggies_to_eat = food_recom.veggies["general_health"]["female"]["eat"]
                veggies_to_avoid = food_recom.veggies["general_health"]["female"]["avoid"]
                carbs_to_eat = food_recom.carbs["general_health"]["female"]["eat"]
                carbs_to_avoid = food_recom.carbs["general_health"]["female"]["avoid"]
                fats_to_eat = food_recom.fats["general_health"]["female"]["eat"]
                fats_to_avoid = food_recom.fats["general_health"]["female"]["avoid"]
                proteins_to_eat = food_recom.proteins["general_health"]["female"]["eat"]
                proteins_to_avoid = food_recom.proteins["general_health"]["female"]["avoid"]
                fruits_to_eat = food_recom.fruits["general_health"]["female"]["eat"]
                fruits_to_avoid = food_recom.fruits["general_health"]["female"]["avoid"]
            else:
                veggies_to_eat = food_recom.veggies["general_health"]["male"]["eat"]
                veggies_to_avoid = food_recom.veggies["general_health"]["male"]["avoid"]
                carbs_to_eat = food_recom.carbs["general_health"]["male"]["eat"]
                carbs_to_avoid = food_recom.carbs["general_health"]["male"]["avoid"]
                fats_to_eat = food_recom.fats["general_health"]["male"]["eat"]
                fats_to_avoid = food_recom.fats["general_health"]["male"]["avoid"]
                proteins_to_eat = food_recom.proteins["general_health"]["male"]["eat"]
                proteins_to_avoid = food_recom.proteins["general_health"]["male"]["avoid"]
                fruits_to_eat = food_recom.fruits["general_health"]["male"]["eat"]
                fruits_to_avoid = food_recom.fruits["general_health"]["male"]["avoid"]
        elif goal == ["QUE-GST-NED-BR"]:
            veggies_to_eat = food_recom.veggies["brain_health"]["eat"]
            veggies_to_avoid = food_recom.veggies["brain_health"]["avoid"]
            carbs_to_eat = food_recom.carbs["brain_health"]["eat"]
            carbs_to_avoid = food_recom.carbs["brain_health"]["avoid"]
            fats_to_eat = food_recom.fats["brain_health"]["eat"]
            fats_to_avoid = food_recom.fats["brain_health"]["avoid"]
            proteins_to_eat = food_recom.proteins["brain_health"]["eat"]
            proteins_to_avoid = food_recom.proteins["brain_health"]["avoid"]
            fruits_to_eat = food_recom.fruits["brain_health"]["eat"]
            fruits_to_avoid = food_recom.fruits["brain_health"]["avoid"]
        elif goal == ["QUE-GST-NED-GF"]:
            veggies_to_eat = food_recom.veggies["gastro_fitness"]["eat"]
            veggies_to_avoid = food_recom.veggies["gastro_fitness"]["avoid"]
            carbs_to_eat = food_recom.carbs["gastro_fitness"]["eat"]
            carbs_to_avoid = food_recom.carbs["gastro_fitness"]["avoid"]
            fats_to_eat = food_recom.fats["gastro_fitness"]["eat"]
            fats_to_avoid = food_recom.fats["gastro_fitness"]["avoid"]
            proteins_to_eat = food_recom.proteins["gastro_fitness"]["eat"]
            proteins_to_avoid = food_recom.proteins["gastro_fitness"]["avoid"]
            fruits_to_eat = food_recom.fruits["gastro_fitness"]["eat"]
            fruits_to_avoid = food_recom.fruits["gastro_fitness"]["avoid"]

        elif goal == ["QUE-GST-NED-PCN"]:
            veggies_to_eat = food_recom.veggies["post_covid"]["eat"]
            veggies_to_avoid = food_recom.veggies["post_covid"]["avoid"]
            carbs_to_eat = food_recom.carbs["post_covid"]["eat"]
            carbs_to_avoid = food_recom.carbs["post_covid"]["avoid"]
            fats_to_eat = food_recom.fats["post_covid"]["eat"]
            fats_to_avoid = food_recom.fats["post_covid"]["avoid"]
            proteins_to_eat = food_recom.proteins["post_covid"]["eat"]
            proteins_to_avoid = food_recom.proteins["post_covid"]["avoid"]
            fruits_to_eat = food_recom.fruits["post_covid"]["eat"]
            fruits_to_avoid = food_recom.fruits["post_covid"]["avoid"]

        elif goal == ["QUE-GST-NED-PE"]:
            veggies_to_eat = food_recom.veggies["performance_and_endurance"]["eat"]
            veggies_to_avoid = food_recom.veggies["performance_and_endurance"]["avoid"]
            carbs_to_eat = food_recom.carbs["performance_and_endurance"]["eat"]
            carbs_to_avoid = food_recom.carbs["performance_and_endurance"]["avoid"]
            fats_to_eat = food_recom.fats["performance_and_endurance"]["eat"]
            fats_to_avoid = food_recom.fats["performance_and_endurance"]["avoid"]
            proteins_to_eat = food_recom.proteins["performance_and_endurance"]["eat"]
            proteins_to_avoid = food_recom.proteins["performance_and_endurance"]["avoid"]
            fruits_to_eat = food_recom.fruits["performance_and_endurance"]["eat"]
            fruits_to_avoid = food_recom.fruits["performance_and_endurance"]["avoid"]

        elif goal == ["QUE-GST-NED-SK"]:
            veggies_to_eat = food_recom.veggies["skin"]["eat"]
            veggies_to_avoid = food_recom.veggies["skin"]["avoid"]
            carbs_to_eat = food_recom.carbs["skin"]["eat"]
            carbs_to_avoid = food_recom.carbs["skin"]["avoid"]
            fats_to_eat = food_recom.fats["skin"]["eat"]
            fats_to_avoid = food_recom.fats["skin"]["avoid"]
            proteins_to_eat = food_recom.proteins["skin"]["eat"]
            proteins_to_avoid = food_recom.proteins["skin"]["avoid"]
            fruits_to_eat = food_recom.fruits["skin"]["eat"]
            fruits_to_avoid = food_recom.fruits["skin"]["avoid"]

        elif goal == ["QUE-GST-NED-SP"]:
            veggies_to_eat = food_recom.veggies["sleep"]["eat"]
            veggies_to_avoid = food_recom.veggies["sleep"]["avoid"]
            carbs_to_eat = food_recom.carbs["sleep"]["eat"]
            carbs_to_avoid = food_recom.carbs["sleep"]["avoid"]
            fats_to_eat = food_recom.fats["sleep"]["eat"]
            fats_to_avoid = food_recom.fats["sleep"]["avoid"]
            proteins_to_eat = food_recom.proteins["sleep"]["eat"]
            proteins_to_avoid = food_recom.proteins["sleep"]["avoid"]
            fruits_to_eat = food_recom.fruits["sleep"]["eat"]
            fruits_to_avoid = food_recom.fruits["sleep"]["avoid"]

        elif goal == ["QUE-GST-NED-SW"]:
            veggies_to_eat = food_recom.veggies["shiftworker"]["eat"]
            veggies_to_avoid = food_recom.veggies["shiftworker"]["avoid"]
            carbs_to_eat = food_recom.carbs["shiftworker"]["eat"]
            carbs_to_avoid = food_recom.carbs["shiftworker"]["avoid"]
            fats_to_eat = food_recom.fats["shiftworker"]["eat"]
            fats_to_avoid = food_recom.fats["shiftworker"]["avoid"]
            proteins_to_eat = food_recom.proteins["shiftworker"]["eat"]
            proteins_to_avoid = food_recom.proteins["shiftworker"]["avoid"]
            fruits_to_eat = food_recom.fruits["shiftworker"]["eat"]
            fruits_to_avoid = food_recom.fruits["shiftworker"]["avoid"]

        if goal == ["QUE-GST-NED-WL"]:
            veggies_to_eat = food_recom.veggies["weight_reduction"]["eat"]
            veggies_to_avoid = food_recom.veggies["weight_reduction"]["avoid"]
            carbs_to_eat = food_recom.carbs["weight_reduction"]["eat"]
            carbs_to_avoid = food_recom.carbs["weight_reduction"]["avoid"]
            fats_to_eat = food_recom.fats["weight_reduction"]["eat"]
            fats_to_avoid = food_recom.fats["weight_reduction"]["avoid"]
            proteins_to_eat = food_recom.proteins["weight_reduction"]["eat"]
            proteins_to_avoid = food_recom.proteins["weight_reduction"]["avoid"]
            fruits_to_eat = food_recom.fruits["weight_reduction"]["eat"]
            fruits_to_avoid = food_recom.fruits["weight_reduction"]["avoid"]

        return veggies_to_eat, veggies_to_avoid, proteins_to_eat, proteins_to_avoid, fats_to_eat, fats_to_avoid, carbs_to_eat, carbs_to_avoid, fruits_to_eat, fruits_to_avoid

    @staticmethod
    def alcohol_intake(intake_data):
        intake_data = round(intake_data, 1)
        if intake_data > 14:
            inference = "*Your intake is beyond the permitted limit."
        elif 0< intake_data <= 14:
            inference = "*Your intake is in the permitted limit."
        elif intake_data == 0 or intake_data < 0 or intake_data == None:
            inference = "*You don't intake any alcohol"
        intake=intake_data    
        return intake, inference

    @staticmethod
    def caffine_intake(intake_data):
        if intake_data > 400:
            inference = "Your caffeine intake is within range*"
        elif intake_data <= 400:
            inference = "Your caffeine intake is within range*"
        elif intake_data == 0 or intake_data < 0 or intake_data == None:
            inference = "Your caffeine intake is 0*"
        intake = round(intake_data, 1)
        return intake, inference

    @staticmethod
    def fluid_intake_inference(fluids):
        max_value = max(zip(fluids.values(), fluids.keys()))[1]
        inference = {
            "tea": "Majority of calorie intake contribution among beverages is from Tea/ Coffee",
            "carbonated": "Majority of calorie intake contribution among beverages is from Carbonated beverages",
            "fruit": "Majority of calorie intake contribution among beverages is from Juices",
            "alc": "Majority of calorie intake contribution among beverages is from alcohol",
            "milk": "Majority of calorie intake contribution among beverages is from Milk",
        }
        return inference[max_value]

    @staticmethod
    def cal_percentage(specific_calorie, total_calorie):
        if total_calorie == 0:
            per = 0
        else:
            per = (specific_calorie / total_calorie) * 100
        return round(per, 1)

    @staticmethod
    def h_graph_calc(nutrition, microbiome=None, genetics=None):
        # microbiome = [51.2,46.64,47.67,46.49,38.54]
        # genetics = [39.29,55,57.69,46.43,37.5]
        nutrition__ = nutrition[:4]
        percentage_list = [22,10,20,48]
        nutrition_list = [nutrition__[i] *100/percentage_list[i] for i in range(len(nutrition__))]
        nutrition_list.append(nutrition[4])
        x_nara_weitage = [.3, .2, .1, .2, .2]
        if microbiome == None and genetics == None:
            weitage_list = nutrition_list
        if microbiome != None and genetics != None and nutrition != None:
            nutrition_list = [i * .7 for i in nutrition_list]
            microbiome_list = [i * .15 for i in microbiome]
            genetics_list = [i * .15 for i in genetics]

            weitage_list = [nutrition_list[i] + microbiome_list[i] + genetics_list[i] for i in
                            range(len(nutrition))]

        if genetics == None and microbiome != None and nutrition != None:
            nutrition_list = [i * .85 for i in nutrition_list]
            microbiome_list = [i * .15 for i in microbiome]

            weitage_list = [nutrition_list[i] + microbiome_list[i]
                            for i in range(len(nutrition))]

        if microbiome == None and genetics != None and nutrition != None:
            nutrition_list = [i * .85 for i in nutrition_list]
            genetics_list = [i * .15 for i in genetics]

            weitage_list = [nutrition_list[i] + genetics_list[i]
                            for i in range(len(nutrition))]

        hgraph_score = sum([weitage_list[i] * x_nara_weitage[i]
                           for i in range(len(weitage_list))])
        return round(hgraph_score,1)

    @staticmethod
    def get_further_recommendations(consumer_internal_markers):
        goal_dict = {
            "QUE-GST-NED-GH":"general_health",
            "QUE-GST-NED-BR":"brain_health",
            "QUE-GST-NED-GF":"gastro_fitness",
            "QUE-GST-NED-PCN":"post_covid",
            "QUE-GST-NED-PE":"performance_and_endurance",
            "QUE-GST-NED-SK":"skin",
            "QUE-GST-NED-SP":"sleep",
            "QUE-GST-NED-WL":"weight_loss",
            "QUE-GST-NED-SW":"shift_worker",
        }
        genetics = []
        microbiome = []
        goal_ = [k for k, v in consumer_internal_markers.items() if v == True]
        goal = goal_dict.get(goal_[0])
        if goal == "post_covid" or goal == "shift_worker":
            return genetics, microbiome
        else:
            genetics = tests_recommend.tests[goal]["genetics"]
            microbiome = tests_recommend.tests[goal]["microbiome"]
            return genetics, microbiome

    @staticmethod
    def personalized_recomm(diabetes, lung_risk, cvd, goal_data):
        goals = {
            "QUE-GST-NED-GH": "General Health",
            "QUE-GST-NED-BR": "Brain Health",
            "QUE-GST-NED-GF": "Gastro Fitness",
            "QUE-GST-NED-PCN": "Post-Covid Nutrition",
            "QUE-GST-NED-PE": "Performance and Endurance",
            "QUE-GST-NED-SK": "Skin",
            "QUE-GST-NED-SP": "Sleep",
            "QUE-GST-NED-SW": "Shift Worker",
            "QUE-GST-NED-WL": "Weight Loss",
        }
        inference = []
        goal = [k for k, v in goal_data.items() if v == True]
        if diabetes > 7:
            inference.append("You have high risk of diabetes")
        elif diabetes < 7:
            inference.append("You have low risk of diabetes")
        if cvd <= 5:
            inference.append(
                "You have very low risk of having cardiovascular diseases")
        elif 5 < cvd <= 10:
            inference.append(
                "You have low risk of having cardiovascular diseases")
        elif 10 < cvd < 30:
            inference.append(
                "You have high risk of having cardiovascular diseases")
        elif cvd >= 30:
            inference.append(
                "You have very high risk of having cardiovascular diseases")
        if lung_risk > 20:
            inference.append("You have very high risk of having Lung diseases")
        elif 0 < lung_risk < 20:
            inference.append("You have risk of having Lung diseases")
        elif lung_risk == 0:
            inference.append("You have very low risk of having Lung diseases")

        return inference, goals[goal[0]]

    @staticmethod
    def get_snack_cal(fru=0, veg=0, nut=0, fried=0, biscut_sweets=0):
        total_cal = fru+veg+nut+fried+biscut_sweets
        return total_cal

    @staticmethod
    def get_insomnia_score(data):
        slpd = {'QUE-GST-SLPD-M': 1, 'QUE-GST-SLPD-MO': 2,
                'QUE-GST-SLPD-N': 0, 'QUE-GST-SLPD-S': 3, 'QUE-GST-SLPD-VS': 4}
        slpsa = {'QUE-GST-SLPSA-M': 1, 'QUE-GST-SLPSA-MO': 2,
                 'QUE-GST-SLPSA-N': 0, 'QUE-GST-SLPSA-S': 3, 'QUE-GST-SLPSA-VS': 4}
        wke = {'QUE-GST-WKE-M': 1, 'QUE-GST-WKE-MO': 2,
               'QUE-GST-WKE-N': 0, 'QUE-GST-WKE-S': 3, 'QUE-GST-WKE-VS': 4}
        slssat = {'QUE-GST-SLSSAT-M': 1, 'QUE-GST-SLSSAT-MO': 2, 'QUE-GST-SLSSAT-N': 0, 'QUE-GST-SLSSAT-S': 3,
                  'QUE-GST-SLSSAT-VS': 4}
        slsnot = {'QUE-GST-SLSNOT-M': 1, 'QUE-GST-SLSNOT-MO': 2, 'QUE-GST-SLSNOT-N': 0, 'QUE-GST-SLSNOT-S': 3,
                  'QUE-GST-SLSNOT-VS': 4}
        slpwor = {'QUE-GST-SLPWOR-M': 1, 'QUE-GST-SLPWOR-MO': 2, 'QUE-GST-SLPWOR-N': 0, 'QUE-GST-SLPWOR-S': 3,
                  'QUE-GST-SLPWOR-VS': 4}
        slpint = {'QUE-GST-SLPINT-M': 1, 'QUE-GST-SLPINT-MO': 2, 'QUE-GST-SLPINT-N': 0, 'QUE-GST-SLPINT-S': 3,
                  'QUE-GST-SLPINT-VS': 4}
        insomnia_score = slpd[data['QUE-GST-SLPD']]+slpsa[data['QUE-GST-SLPSA']]+wke[data['QUE-GST-WKE']] +\
            slssat[data['QUE-GST-SLSSAT']]+slsnot[data['QUE-GST-SLSNOT']]+slpwor[data['QUE-GST-SLPWOR']] +\
            slpint[data['QUE-GST-SLPINT']]

        return insomnia_score

    @staticmethod
    def get_osa_score(data):
        apnea_score = 0
        gender = 'F' if data.get("QUE-GIN-GEN-F") else 'M'
        weight = data.get("QUE-GST-CWG-KG", 0)
        height = data.get("QUE-GST-HGT-CM", 0)
        bmi = Calculate.calc_bmi(height, weight)

        const_dict = {
            'QUE-GST-SNRLOU-Y-': True,
            "QUE-GST-SLPTIR-Y-": True,
            "QUE-GST-SLPSTP-Y-": True,
            "QUE-GST-SLPBLD-Y-": True,
            "QUE-GST-SLPNCK-Y-": True,
        }
        all_que_list = ["QUE-GST-SNRLOU", "QUE-GST-SLPTIR",
                        "QUE-GST-SLPSTP", "QUE-GST-SLPBLD", "QUE-GST-SLPNCK"]
        for i in all_que_list:
            if const_dict.get(data[i], False):
                apnea_score += 1

        if gender == "M":
            apnea_score += 1
        if bmi > 35:
            apnea_score += 1
        if data.get("QUE-GIN-AGE") > 50:
            apnea_score += 1

        return apnea_score

    @staticmethod
    def get_rls_score(data):
        rlsams = {'QUE-GST-RLSAMS-M': 2, 'QUE-GST-RLSAMS-MD': 1, 'QUE-GST-RLSAMS-N-': 0, 'QUE-GST-RLSAMS-S': 3,
                  'QUE-GST-RLSAMS-V': 4}
        rlsams_mov = {'QUE-GST-RLSAMS-MOV-M': 2, 'QUE-GST-RLSAMS-MOV-MD': 1, 'QUE-GST-RLSAMS-MOV-N-': 0,
                      'QUE-GST-RLSAMS-MOV-S': 3, 'QUE-GST-RLSAMS-MOV-V': 4}
        rlsds = {'QUE-GST-RLSDS-MOV-M': 2, 'QUE-GST-RLSDS-MOV-MD': 1, 'QUE-GST-RLSDS-MOV-N-': 0,
                 'QUE-GST-RLSDS-MOV-S': 3, 'QUE-GST-RLSDS-MOV-V': 4}
        rlsdsho_mov = {'QUE-GST-RLSDSHO-MOV-M': 2, 'QUE-GST-RLSDSHO-MOV-MD': 1, 'QUE-GST-RLSDSHO-MOV-N-': 0,
                       'QUE-GST-RLSDSHO-MOV-S': 3, 'QUE-GST-RLSDSHO-MOV-V': 4}
        rlsdsho_ov = {'QUE-GST-RLSDSHO-OV-M': 2, 'QUE-GST-RLSDSHO-OV-MD': 1, 'QUE-GST-RLSDSHO-OV-N-': 0,
                      'QUE-GST-RLSDSHO-OV-S': 3, 'QUE-GST-RLSDSHO-OV-V': 4}
        rlsdsho_rls = {'QUE-GST-RLSDSHO-RLS-M': 2, 'QUE-GST-RLSDSHO-RLS-MD': 1, 'QUE-GST-RLSDSHO-RLS-N-': 0,
                       'QUE-GST-RLSDSHO-RLS-V': 4, 'QUE-GST-RLSDSHO-RLS-S': 3}
        rlsdsho_rls_sev = {'QUE-GST-RLSDSHO-RLS-SEV-M': 2, 'QUE-GST-RLSDSHO-RLS-SEV-N-': 0,
                           'QUE-GST-RLSDSHO-RLS-SEV-S': 3, 'QUE-GST-RLSDSHO-RLS-SEV-V': 4,
                           'QUE-GST-RLSDSHO-RLS-SEV-MD': 1}
        rlsdsho_rls_sev_imp = {'QUE-GST-RLSDSHO-RLS-SEV-IMP-M': 2, 'QUE-GST-RLSDSHO-RLS-SEV-IMP-MD': 1,
                               'QUE-GST-RLSDSHO-RLS-SEV-IMP-N-': 0, 'QUE-GST-RLSDSHO-RLS-SEV-IMP-S': 3,
                               'QUE-GST-RLSDSHO-RLS-SEV-IMP-V': 4}
        rlsdsho_rls_sev_mda = {'QUE-GST-RLSDSHO-RLS-SEV-MDA-M': 2, 'QUE-GST-RLSDSHO-RLS-SEV-MDA-MD': 1,
                               'QUE-GST-RLSDSHO-RLS-SEV-MDA-N-': 0, 'QUE-GST-RLSDSHO-RLS-SEV-MDA-S': 3,
                               'QUE-GST-RLSDSHO-RLS-SEV-MDA-V': 4}
        rlsdsho_rls_sev_rls = {'QUE-GST-RLSDSHO-RLS-SEV-RLS-M': 2, 'QUE-GST-RLSDSHO-RLS-SEV-RLS-MD': 1,
                               'QUE-GST-RLSDSHO-RLS-SEV-RLS-N-': 0, 'QUE-GST-RLSDSHO-RLS-SEV-RLS-S': 3,
                               'QUE-GST-RLSDSHO-RLS-SEV-RLS-V': 4}

        rls_score = sum([
            rlsams[data["QUE-GST-RLSAMS"]],
            rlsams_mov[data['QUE-GST-RLSAMS-MOV']],
            rlsds[data['QUE-GST-RLSDS']],
            rlsdsho_mov[data['QUE-GST-RLSDSHO']],
            rlsdsho_ov[data['QUE-GST-RLSDSHO-OV']],
            rlsdsho_rls[data["QUE-GST-RLSDSHO-RLS"]],
            rlsdsho_rls_sev[data['QUE-GST-RLSDSHO-RLS-SEV']],
            rlsdsho_rls_sev_imp[data['QUE-GST-RLSDSHO-RLS-SEV-IMP']],
            rlsdsho_rls_sev_mda[data['QUE-GST-RLSDSHO-RLS-SEV-MDA']],
            rlsdsho_rls_sev_rls[data["QUE-GST-RLSDSHO-RLS-SEV-RLS"]]
        ]
        )

        return rls_score

    @staticmethod
    def get_stress_score(data):
        cntr = {'QUE-GST-STR-CNTR-NV': 0, 'QUE-GST-STR-CNTR-ANV': 1, 'QUE-GST-STR-CNTR-ST': 2, 'QUE-GST-STR-CNTR-FO': 3,
                'QUE-GST-STR-CNTR-VO': 4}
        conf = {'QUE-GST-STR-CONF-NV': 4, 'QUE-GST-STR-CONF-ANV': 3, 'QUE-GST-STR-CONF-ST': 2, 'QUE-GST-STR-CONF-FO': 1,
                'QUE-GST-STR-CONF-VO': 0}
        yway = {'QUE-GST-STR-YWAY-NV': 4, 'QUE-GST-STR-YWAY-ANV': 3, 'QUE-GST-STR-YWAY-ST': 2, 'QUE-GST-STR-YWAY-FO': 1,
                'QUE-GST-STR-YWAY-VO': 0}
        diff = {'QUE-GST-STR-DIFF-NV': 0, 'QUE-GST-STR-DIFF-ANV': 1, 'QUE-GST-STR-DIFF-ST': 2, 'QUE-GST-STR-DIFF-FO': 3,
                'QUE-GST-STR-DIFF-VO': 4}

        stress_score = sum([cntr[data["QUE-GST-STR-CNTR"]], conf[data['QUE-GST-STR-CONF']],
                            yway[data['QUE-GST-STR-YWAY']], diff[data['QUE-GST-STR-DIFF']]])
        return stress_score

    @staticmethod
    def get_anxiety_score(data):
        afr = {'QUE-GST-PR-PRB-PRB-AFR-PSO': 0, 'QUE-GST-PR-PRB-PRB-AFR-DER': 1, 'QUE-GST-PR-PRB-PRB-AFR-VIT': 2,
               'QUE-GST-PR-PRB-PRB-AFR-VITN': 3}
        con = {'QUE-GST-PR-PRB-PRB-CON-PSO': 0, 'QUE-GST-PR-PRB-PRB-CON-DER': 1, 'QUE-GST-PR-PRB-PRB-CON-VIT': 2,
               'QUE-GST-PR-PRB-PRB-CON-VITN': 3}
        irr = {'QUE-GST-PR-PRB-PRB-IRR-PSO': 0, 'QUE-GST-PR-PRB-PRB-IRR-DER': 1, 'QUE-GST-PR-PRB-PRB-IRR-VIT': 2,
               'QUE-GST-PR-PRB-PRB-IRR-VITN': 3}
        nar = {'QUE-GST-PR-PRB-PRB-NAR-PSO': 0, 'QUE-GST-PR-PRB-PRB-NAR-DER': 1, 'QUE-GST-PR-PRB-PRB-NAR-VIT': 2,
               'QUE-GST-PR-PRB-PRB-NAR-VITN': 3}
        rel = {'QUE-GST-PR-PRB-PRB-REL-PSO': 0, 'QUE-GST-PR-PRB-PRB-REL-DER': 1, 'QUE-GST-PR-PRB-PRB-REL-VIT': 2,
               'QUE-GST-PR-PRB-PRB-REL-VITN': 3}
        rst = {'QUE-GST-PR-PRB-PRB-RST-PSO': 0, ' QUE-GST-PR-PRB-PRB-RST-DER': 1, 'QUE-GST-PR-PRB-PRB-RST-VIT': 2,
               'QUE-GST-PR-PRB-PRB-RST-VITN': 3}
        woo = {'QUE-GST-PR-PRB-PRB-WOO-PSO': 0, 'QUE-GST-PR-PRB-PRB-WOO-DER': 1, 'QUE-GST-PR-PRB-PRB-WOO-VIT': 2,
               'QUE-GST-PR-PRB-PRB-WOO-VITN': 3}

        anxiety_score = sum([
            afr[data["QUE-GST-PR-PRB-PRB-AFR"]],
            con[data["QUE-GST-PR-PRB-PRB-CON"]],
            irr[data["QUE-GST-PR-PRB-PRB-IRR"]],
            nar[data["QUE-GST-PR-PRB-PRB-NAR"]],
            rel[data["QUE-GST-PR-PRB-PRB-REL"]],
            rst[data["QUE-GST-PR-PRB-PRB-RST"]],
            woo[data["QUE-GST-PR-PRB-PRB-WOO"]]
        ]
        )

        return anxiety_score

    @staticmethod
    def get_depression_score(data):
        apt = {'QUE-GST-PR-PRB-APT-PSO': 0, 'QUE-GST-PR-PRB-APT-DER': 1, 'QUE-GST-PR-PRB-APT-VIT': 2,
               'QUE-GST-PR-PRB-APT-VITN': 3}
        fed = {'QUE-GST-PR-PRB-FED-PSO': 0, 'QUE-GST-PR-PRB-FED-DER': 1, 'QUE-GST-PR-PRB-FED-VIT': 2,
               'QUE-GST-PR-PRB-FED-VITN': 3}
        hur = {'QUE-GST-PR-PRB-HUR-PSO': 0, 'QUE-GST-PR-PRB-HUR-DER': 1, 'QUE-GST-PR-PRB-HUR-VIT': 2,
               'QUE-GST-PR-PRB-HUR-VITN': 3}
        lit = {'QUE-GST-PR-PRB-LIT-PSO': 0, 'QUE-GST-PR-PRB-LIT-DER': 1, 'QUE-GST-PR-PRB-LIT-VIT': 2,
               'QUE-GST-PR-PRB-LIT-VITN': 3}
        mov = {'QUE-GST-PR-PRB-MOV-PSO': 0, 'QUE-GST-PR-PRB-MOV-DER': 1, 'QUE-GST-PR-PRB-MOV-VIT': 2,
               'QUE-GST-PR-PRB-MOV-VITN': 3}
        sel = {'QUE-GST-PR-PRB-SEL-PSO': 0, 'QUE-GST-PR-PRB-SEL-DER': 1, 'QUE-GST-PR-PRB-SEL-VIT': 2,
               'QUE-GST-PR-PRB-SEL-VITN': 3}
        tas = {'QUE-GST-PR-PRB-TAS-PSO': 0, 'QUE-GST-PR-PRB-TAS-DER': 1, 'QUE-GST-PR-PRB-TAS-VIT': 2,
               'QUE-GST-PR-PRB-TAS-VITN': 3}
        toc = {'QUE-GST-PR-PRB-TOC-PSO': 0, 'QUE-GST-PR-PRB-TOC-DER': 1, 'QUE-GST-PR-PRB-TOC-VIT': 2,
               'QUE-GST-PR-PRB-TOC-VITN': 3}
        tre = {'QUE-GST-PR-PRB-TRE-PSO': 0, 'QUE-GST-PR-PRB-TRE-DER': 1, 'QUE-GST-PR-PRB-TRE-VIT': 2,
               'QUE-GST-PR-PRB-TRE-VITN': 3}

        depression_score = sum([
            apt[data["QUE-GST-PR-PRB-APT"]],
            fed[data["QUE-GST-PR-PRB-FED"]],
            hur[data["QUE-GST-PR-PRB-HUR"]],
            lit[data["QUE-GST-PR-PRB-LIT"]],
            mov[data["QUE-GST-PR-PRB-MOV"]],
            sel[data["QUE-GST-PR-PRB-SEL"]],
            tas[data["QUE-GST-PR-PRB-TAS"]],
            toc[data["QUE-GST-PR-PRB-TOC"]],
            tre[data["QUE-GST-PR-PRB-TRE"]]
        ]
        )

        return depression_score

    @staticmethod
    def get_nutrition_score(goal, consumer_internal_markers, ethnicity, sbp, cholesterol, spe, sleep_duration, pal, bmi, caffine, fat, protein, grains, alcohol, veg_and_fru_servings):
        nutrition_score = []
        goal = [k for k, v in goal.items() if v == True]
        age = consumer_internal_markers.get("QUE-GIN-AGE", 0)
        gender = 'F' if consumer_internal_markers.get("QUE-GIN-GEN-F") else 'M'
        ethnicity = [k for k, v in ethnicity.items() if v == True]
        sbp = [k for k, v in sbp.items() if v == True]
        cholesterol = [k for k, v in cholesterol.items() if v == True]
        diabetic = consumer_internal_markers.get("QUE-GST-CON-YD")
        physiology_score = 0  # 3pts max physiology
        if 18 <= age <= 30:
            physiology_score += 3
        elif 31 <= age <= 50:
            physiology_score += 2
        elif 60 <= age <= 74:
            physiology_score += 1
        elif age > 75:
            physiology_score += 0

        if diabetic:  # 5 pts max  physiology
            physiology_score += 0
        else:
            physiology_score += 5

        if ethnicity == ["QUE-GIN-ETH-AFR"]:  # bmi classification max 10 pts physiology
            if bmi < 21.9:
                physiology_score += 5
            if 22 <= bmi <= 27.9:
                physiology_score += 10
            if 28 <= bmi <= 32.9:
                physiology_score += 5
            if bmi > 33:
                physiology_score += 0

        elif ethnicity == ["QUE-GIN-ETH-CAU"]:
            if bmi < 18:
                physiology_score += 5
            if 18 <= bmi <= 24.9:
                physiology_score += 10
            if 25 <= bmi <= 29.9:
                physiology_score += 5
            if bmi > 30:
                physiology_score += 0

        elif ethnicity == ["QUE-GIN-ETH-EAS"]:
            if bmi < 18.5:
                physiology_score += 5
            if 18.5 <= bmi <= 23:
                physiology_score += 10
            if 23 <= bmi <= 25:
                physiology_score += 5
            if bmi >= 25:
                physiology_score += 0

        elif ethnicity == ["QUE-GIN-ETH-HIS"] or ethnicity == ["QUE-GIN-ETH-NTA"] or ethnicity == ["QUE-GIN-ETH-PAI"]:
            if bmi < 18.4:
                physiology_score += 5
            if 18.5 <= bmi <= 24.9:
                physiology_score += 10
            if 25 <= bmi <= 29.9:
                physiology_score += 5
            if bmi >= 30:
                physiology_score += 0

        elif ethnicity == ["QUE-GIN-ETH-MID"]:
            if bmi < 18:
                physiology_score += 5
            if 18 <= bmi <= 22.9:
                physiology_score += 10
            if 23 <= bmi <= 27.5:
                physiology_score += 5
            if bmi >= 27.6:
                physiology_score += 0

        elif ethnicity == ["QUE-GIN-ETH-SEA"]:
            if bmi < 18:
                physiology_score += 5
            if 18 <= bmi <= 24.9:
                physiology_score += 10
            if 25 <= bmi <= 29.9:
                physiology_score += 5
            if bmi >= 30:
                physiology_score += 0

        # systolic BP 1 pts max physiology
        if sbp == ["QUE-GST-SBP-BET-140-159"] or sbp == ["QUE-GST-SBP-BET-160-179"] or sbp == ["QUE-GST-SBP-BET-160-179"]:
            physiology_score += 0
        elif sbp == ["QUE-GST-SBP-BET-120-139"]:
            physiology_score += 1
        elif sbp == ["QUE-GST-SBP-BET-90-120"]:
            physiology_score += 1
        elif sbp == ["QUE-GST-SBP-LT-90"]:
            physiology_score += 0
        else:
            physiology_score += 0

        # cholesterol 3 pts max physiology
        if cholesterol == ["QUE-GST-CHL-LT4"]:
            physiology_score += 3
        elif cholesterol == ["QUE-GST-CHL-BET-4-4.9"]:
            physiology_score += 3
        elif cholesterol == ["QUE-GST-CHL-BET-5-5.9"]:
            physiology_score += 2
        elif cholesterol == ["QUE-GST-CHL-BET-6-6.9"]:
            physiology_score += 1
        elif cholesterol == ["QUE-GST-CHL-GT7"]:
            physiology_score += 0
        else:
            physiology_score += 0

        # smoke 3 points max overall health
        overall_score = 0
        if consumer_internal_markers.get("QUE-GST-SMK-Y"):
            overall_score += 0
        elif consumer_internal_markers.get("QUE-GST-SMK-X"):
            overall_score += 1
        elif consumer_internal_markers.get("QUE-GST-SMK-N"):
            overall_score += 3

        # grains intaker percentage  max pts 5
        
        grain_list = ["QUE-HEH-MEA-BF-FP-GR","QUE-HEH-MEA-LN-FP-GR","QUE-HEH-MEA-DN-FP-GR"]
        grain_percentage = 0
        for item in grain_list:
            if consumer_internal_markers.get(item) != None:
                grain_percentage += consumer_internal_markers.get(item)
        grains = grain_percentage / len(grain_list)
        if grains > 20.2:
            overall_score += 5
        elif grains < 20:
            overall_score += 0

        # fat intake percentage max score 10 pts
        if 25 < fat < 35:
            overall_score += 10
        elif fat < 25:
            overall_score += 5
        elif fat > 35:
            overall_score += 0

        # overll health alcohol consumption (daily intake/8 = drinks)  max 10 pts
        if alcohol != None:
            alcohol = alcohol/8
            if alcohol == 0:
                overall_score += 10
            elif gender == "M":
                if 1 <= alcohol <= 2:
                    overall_score += 5
                elif alcohol > 2:
                    overall_score += 0
            elif gender == "F":
                if 0 < alcohol <= 1:
                    overall_score += 5
                elif alcohol > 1:
                    overall_score += 0

        #    protein overallhealth   max pts 10
        protein_list = ["QUE-HEH-MEA-BF-FP-PR","QUE-HEH-MEA-LN-FP-PR","QUE-HEH-MEA-DN-FP-PR"]
        protein_percentage = 0
        for item in protein_list:
            if consumer_internal_markers.get(item) != None:
                protein_percentage += consumer_internal_markers.get(item)
        protein = protein_percentage / len(protein_list)
        if 10 < protein < 35:
            overall_score += 10
        elif 35 > protein < 40:
            overall_score += 5
        elif protein > 40:
            overall_score += 0
        else:
            overall_score += 0
        # water_consumption 0verallhealth max 5 pts
        water_intake = consumer_internal_markers.get("QUE-HEH-LIQ-WA-QN", 0)
        if water_intake != None:
            water_intake = int(water_intake)
            if gender == "F":
                if water_intake != 0:
                    if water_intake > 2700:
                        overall_score += 5
                    elif 2000 <= water_intake < 2700:
                        overall_score += 3
                    elif water_intake < 2000:
                        overall_score += 0
            elif gender == "M":
                if water_intake != 0:
                    if water_intake >= 3700:
                        overall_score += 5
                    elif 2500 <= water_intake < 3700:
                        overall_score += 3
                    elif water_intake < 2500:
                        overall_score += 0

        # (vegetable and fruit intake/100 = servings) max pts 10
        if veg_and_fru_servings > 5:
            overall_score += 10
        elif 2 <= veg_and_fru_servings < 5:
            overall_score += 5
        elif veg_and_fru_servings < 2:
            overall_score += 0

        aesthetic_score = 0
        # Aesthetic PAL value max score 20
        if pal < 1.20:
            aesthetic_score += 0
        elif 1.20 < pal < 1.375:
            aesthetic_score += 5
        elif 1.375 < pal < 1.55:
            aesthetic_score += 10
        elif 1.55 < pal < 1.725:
            aesthetic_score += 15
        elif pal > 1.725:
            aesthetic_score += 20

        # Mental health
        mental_health_score = 0
        # sleep efficiency max score 5pts
        if spe >= 85:
            mental_health_score += 5
        elif 70 <= spe <= 84.99:
            mental_health_score += 3
        elif spe < 70:
            mental_health_score += 0
        # sleep_duration max pts 5 mental health
        if sleep_duration < 7:
            mental_health_score += 0
        elif 7 <= sleep_duration <= 9:
            mental_health_score += 5
        elif sleep_duration > 9:
            mental_health_score += 0
        # Q-17
        nutrition_score.append(physiology_score)
        nutrition_score.append(mental_health_score)
        nutrition_score.append(aesthetic_score)
        nutrition_score.append(overall_score)

    # Goal risk scoring
        if goal == ["QUE-GST-NED-GH"]:
            general_health = physiology_score + \
                mental_health_score + aesthetic_score + overall_score
            nutrition_score.append(general_health)

        if goal == ["QUE-GST-NED-BR"]:
            mental_score_A = 0
            # question number 17# max 10 points
            if consumer_internal_markers.get("QUE-GST-FL-RES-ORG-N"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-FL-RES-DIS-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-FL-RES-ATT-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-FL-RES-WT-N"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-FL-RES-LOST-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-FL-RES-PLS-N"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-FL-RES-GT-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-FL-RES-STM-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-FL-RES-MEM-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-FL-RES-RDB-N"):
                mental_score_A += 1

            # option 2 selected question no 18 max 13 pts
            if consumer_internal_markers.get("QUE-GST-SL-RES-AGO-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SL-RES-ENL-N"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SL-RES-HP-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SL-RES-NGT-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SL-RES-NOT-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SL-RES-PDW-N"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SL-RES-SLF-N"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SL-RES-STC-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SL-RES-STG-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SL-RES-STW-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SL-RES-TEJ-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SL-RES-THS-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SL-RES-WNG-Y"):
                mental_score_A += 1

    # supplimentary set questions
            mental_score_B = 0
            if mental_score_A >= 3:
                important_things_dict = {
                    "QUE-GST-STR-CNTR-NV": 0,
                    "QUE-GST-STR-CNTR-ANV": 1,
                    "QUE-GST-STR-CNTR-ST": 2,
                    "QUE-GST-STR-CNTR-FO": 3,
                    "QUE-GST-STR-CNTR-VO": 4,
                }
                if consumer_internal_markers.get("QUE-GST-STR-CNTR") != None:
                    mental_score_B += important_things_dict.get(
                        consumer_internal_markers.get("QUE-GST-STR-CNTR"), 0)

                confident_dict = {
                    "QUE-GST-STR-CONF-NV": 4,
                    "QUE-GST-STR-CONF-ANV": 3,
                    "QUE-GST-STR-CONF-ST": 2,
                    "QUE-GST-STR-CONF-FO": 1,
                    "QUE-GST-STR-CONF-VO": 0,
                }
                # confident
                if consumer_internal_markers.get("QUE-GST-STR-CONF") != None:
                    mental_score_B += confident_dict.get(
                        consumer_internal_markers.get("QUE-GST-STR-CONF"), 0)

                going_your_way_dict = {
                    "QUE-GST-STR-YWAY-NV": 4,
                    "QUE-GST-STR-YWAY-ANV": 3,
                    "QUE-GST-STR-YWAY-ST": 2,
                    "QUE-GST-STR-YWAY-FO": 1,
                    "QUE-GST-STR-YWAY-VO": 0,
                }
                if consumer_internal_markers.get("QUE-GST-STR-YWAY") != None:
                    mental_score_B += going_your_way_dict.get(
                        consumer_internal_markers.get("QUE-GST-STR-YWAY"), 0)

                difficulties_dict = {
                    "QUE-GST-STR-DIFF-NV": 0,
                    "QUE-GST-STR-DIFF-ANV": 1,
                    "QUE-GST-STR-DIFF-ST": 2,
                    "QUE-GST-STR-DIFF-FO": 3,
                    "QUE-GST-STR-DIFF-VO": 4,
                }
                if consumer_internal_markers.get("QUE-GST-STR-DIFF") != None:
                    mental_score_B += difficulties_dict.get(
                        consumer_internal_markers.get("QUE-GST-STR-DIFF"), 0)
                stress_score = mental_score_B

                # anxienty
                anxiety_1_dict = {
                    "QUE-GST-PR-PRB-PRB-AFR-PSO": 0,
                    "QUE-GST-PR-PRB-PRB-AFR-DER": 1,
                    "QUE-GST-PR-PRB-PRB-AFR-VIT": 2,
                    "QUE-GST-PR-PRB-PRB-AFR-VITN": 3

                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-PRB-AFR") != None:
                    mental_score_B += anxiety_1_dict.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-PRB-AFR"), 0)

                anxiety_2_dict = {
                    "QUE-GST-PR-PRB-PRB-CON-PSO": 0,
                    "QUE-GST-PR-PRB-PRB-CON-DER": 1,
                    "QUE-GST-PR-PRB-PRB-CON-VIT": 2,
                    "QUE-GST-PR-PRB-PRB-CON-VITN": 3,
                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-PRB-CON") != None:
                    mental_score_B += anxiety_2_dict.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-PRB-CON"), 0)

                anxiety_3_dict = {
                    "QUE-GST-PR-PRB-PRB-IRR-PSO": 0,
                    "QUE-GST-PR-PRB-PRB-IRR-DER": 1,
                    "QUE-GST-PR-PRB-PRB-IRR-VIT": 2,
                    "QUE-GST-PR-PRB-PRB-IRR-VITN": 3,
                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-PRB-IRR") != None:
                    mental_score_B += anxiety_3_dict.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-PRB-IRR"), 0)

                anxiety_4_dict = {
                    "QUE-GST-PR-PRB-PRB-NAR-PSO": 0,
                    "QUE-GST-PR-PRB-PRB-NAR-DER": 1,
                    "QUE-GST-PR-PRB-PRB-NAR-VIT": 2,
                    "QUE-GST-PR-PRB-PRB-NAR-VITN": 3,
                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-PRB-NAR") != None:
                    mental_score_B += anxiety_4_dict.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-PRB-NAR"), 0)

                anxiety_5_dict = {
                    "QUE-GST-PR-PRB-PRB-REL-PSO": 0,
                    "QUE-GST-PR-PRB-PRB-REL-DER": 1,
                    "QUE-GST-PR-PRB-PRB-REL-VIT": 2,
                    "QUE-GST-PR-PRB-PRB-REL-VITN": 3,
                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-PRB-REL") != None:
                    mental_score_B += anxiety_5_dict.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-PRB-REL"), 0)

                anxiety_6_dict = {
                    "QUE-GST-PR-PRB-PRB-RST-PSO": 0,
                    "QUE-GST-PR-PRB-PRB-RST-DER": 1,
                    "QUE-GST-PR-PRB-PRB-RST-VIT": 2,
                    "QUE-GST-PR-PRB-PRB-RST-VITN": 3,
                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-PRB-RST") != None:
                    mental_score_B += anxiety_6_dict.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-PRB-RST"), 0)

                anxiety_7_dict = {
                    "QUE-GST-PR-PRB-PRB-WOO-PSO": 0,
                    "QUE-GST-PR-PRB-PRB-WOO-DER": 1,
                    "QUE-GST-PR-PRB-PRB-WOO-VIT": 2,
                    "QUE-GST-PR-PRB-PRB-WOO-VITN": 3,
                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-PRB-WOO") != None:
                    mental_score_B += anxiety_7_dict.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-PRB-WOO"), 0)

                # depression scale
                depression_scale_1 = {
                    "QUE-GST-PR-PRB-APT-PSO": 0,
                    "QUE-GST-PR-PRB-APT-DER": 1,
                    "QUE-GST-PR-PRB-APT-VIT": 2,
                    "QUE-GST-PR-PRB-APT-VITN": 3,
                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-APT") != None:
                    mental_score_B += depression_scale_1.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-APT"), 0)

                depression_scale_2 = {
                    "QUE-GST-PR-PRB-FED-PSO": 0,
                    "QUE-GST-PR-PRB-FED-DER": 1,
                    "QUE-GST-PR-PRB-FED-VIT": 2,
                    "QUE-GST-PR-PRB-FED-VITN": 3,
                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-FED") != None:
                    mental_score_B += depression_scale_2.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-FED"), 0)

                depression_scale_3 = {
                    "QUE-GST-PR-PRB-HUR-PSO": 0,
                    "QUE-GST-PR-PRB-HUR-DER": 1,
                    "QUE-GST-PR-PRB-HUR-VIT": 2,
                    "QUE-GST-PR-PRB-HUR-VITN": 3,
                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-HUR") != None:
                    mental_score_B += depression_scale_3.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-HUR"), 0)

                depression_scale_4 = {
                    "QUE-GST-PR-PRB-LIT-PSO": 0,
                    "QUE-GST-PR-PRB-LIT-DER": 1,
                    "QUE-GST-PR-PRB-LIT-VIT": 2,
                    "QUE-GST-PR-PRB-LIT-VITN": 3,
                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-LIT") != None:
                    mental_score_B += depression_scale_4.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-LIT"), 0)

                depression_scale_5 = {
                    "QUE-GST-PR-PRB-MOV-PSO": 0,
                    "QUE-GST-PR-PRB-MOV-DER": 1,
                    "QUE-GST-PR-PRB-MOV-VIT": 2,
                    "QUE-GST-PR-PRB-MOV-VITN": 3,
                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-MOV") != None:
                    mental_score_B += depression_scale_5.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-MOV"), 0)

                depression_scale_6 = {
                    "QUE-GST-PR-PRB-SEL-PSO": 0,
                    "QUE-GST-PR-PRB-SEL-DER": 1,
                    "QUE-GST-PR-PRB-SEL-VIT": 2,
                    "QUE-GST-PR-PRB-SEL-VITN": 3,
                }

                depression_scale_7 = {
                    "QUE-GST-PR-PRB-TAS-PSO": 0,
                    "QUE-GST-PR-PRB-TAS-DER": 1,
                    "QUE-GST-PR-PRB-TAS-VIT": 2,
                    "QUE-GST-PR-PRB-TAS-VITN": 3,
                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-TAS") != None:
                    mental_score_B += depression_scale_7.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-TAS"), 0)

                depression_scale_8 = {
                    "QUE-GST-PR-PRB-TOC-PSO": 0,
                    "QUE-GST-PR-PRB-TOC-DER": 1,
                    "QUE-GST-PR-PRB-TOC-VIT": 2,
                    "QUE-GST-PR-PRB-TOC-VITN": 3,
                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-TOC") != None:
                    mental_score_B += depression_scale_8.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-TOC"), 0)

                depression_scale_9 = {
                    "QUE-GST-PR-PRB-TRE-PSO": 0,
                    "QUE-GST-PR-PRB-TRE-DER": 1,
                    "QUE-GST-PR-PRB-TRE-VIT": 2,
                    "QUE-GST-PR-PRB-TRE-VITN": 3,
                }
                if consumer_internal_markers.get("QUE-GST-PR-PRB-TRE") != None:
                    mental_score_B += depression_scale_9.get(
                        consumer_internal_markers.get("QUE-GST-PR-PRB-TRE"), 0)

            #   18.1   sleep related complaints max 6 points
            if consumer_internal_markers.get("QUE-GST-SLP-CMP-INS-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SLP-CMP-SBR-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SLP-CMP-SLPD-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SLP-CMP-SNR-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SLP-CMP-SOD-Y"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SLP-CMP-TFS-Y"):
                mental_score_A += 1

            # set c move to under the supplimentry questions
            very_vigrous = 0 if consumer_internal_markers.get(
                "QUE-HEH-EAL-EA-VVE-M") == None else consumer_internal_markers.get("QUE-HEH-EAL-EA-VVE-M", 0)
            vigrous = 0 if consumer_internal_markers.get(
                "QUE-HEH-EAL-VA-VE-M") == None else consumer_internal_markers.get("QUE-HEH-EAL-VA-VE-M", 0)
            light = 0 if consumer_internal_markers.get(
                "QUE-HEH-EAL-LG-LE-M") == None else consumer_internal_markers.get("QUE-HEH-EAL-LG-LE-M", 0)
            moderate = 0 if consumer_internal_markers.get(
                "QUE-HEH-EAL-MD-ME-M") == None else consumer_internal_markers.get("QUE-HEH-EAL-MD-ME-M", 0)
            if very_vigrous + vigrous + light + moderate <= 150:
                mental_score_A += 1

            if sleep_duration <= 7:
                mental_score_A += 1

            tea_coffee = 0 if consumer_internal_markers.get(
                "QUE-HEH-LIQ-CP-TC") == None else consumer_internal_markers.get("QUE-HEH-LIQ-CP-TC", 0)
            if tea_coffee * 100/30 >= 8:
                mental_score_A += 1

            # sweetend drinks everyday
            if consumer_internal_markers.get("QUE-HEH-LIQ-SC"):
                mental_score_A += 1

            if consumer_internal_markers.get("QUE-GST-SMK-Y"):
                mental_score_A += 1
            try:
                if stress_score > 8:
                    mental_score_A += 1
            except:
                mental_score_A += 0

            beer_qu = 0 if consumer_internal_markers.get(
                "QUE-HEH-LIQ-ALC-BE") == None else consumer_internal_markers.get("QUE-HEH-LIQ-ALC-BE")/194.117647059
            wine_qu = 0 if consumer_internal_markers.get(
                "QUE-HEH-LIQ-ALC-WN") == None else consumer_internal_markers.get("QUE-HEH-LIQ-ALC-WN")/83.333333333
            liqr_qu = 0 if consumer_internal_markers.get(
                "QUE-HEH-LIQ-ALC-HL") == None else consumer_internal_markers.get("QUE-HEH-LIQ-ALC-HL")/25

            alc_unit = beer_qu + wine_qu + liqr_qu
            if alc_unit >= 14:
                mental_score_A += 1

            try:

                total_mental_score = (mental_score_A + mental_score_B)/(39+64)
            except:
                total_mental_score = mental_score_A/39

            risk_score_percentage = total_mental_score*100
            mental_health_score = 100 - risk_score_percentage
            nutrition_score.append(mental_health_score)

        if goal == ["QUE-GST-NED-PE"]:

            endurence_score = 0
            if consumer_internal_markers.get("QUE-GST-ATH-TY-BG"):
                endurence_score += 40
            if consumer_internal_markers.get("QUE-GST-ATH-TY-TR"):
                endurence_score += 30
            if consumer_internal_markers.get("QUE-GST-ATH-TY-ATR"):
                endurence_score += 20
            if consumer_internal_markers.get("QUE-GST-ATH-TY-PRO"):
                endurence_score += 10

            endurence_score = (endurence_score/40)/2

            pe_score = 0
            ws = int(consumer_internal_markers.get("QUE-GST-HRS-WS-CM"))

            if gender == 'F':
                if ws < 80:
                    pe_score += 0
                elif 80 < ws < 90:
                    pe_score += 5
                elif 90 < ws < 100:
                    pe_score += 100
                elif ws >= 100:
                    pe_score += 15

            if gender == 'M':
                if ws < 90:
                    pe_score += 0
                if 90 <= ws <= 100:
                    pe_score += 5
                if 100 < ws <= 110:
                    pe_score += 10
                if 100 < ws <= 120:
                    pe_score += 15

            if bmi > 25:
                pe_score += 5
            if bmi > 30:
                pe_score += 10
            if 20 < bmi < 25:
                pe_score += 0
            if 18.5 < bmi < 20:
                pe_score += 5
            if bmi < 18.5:
                pe_score += 10

            body_fat_dict = {
                "QUE-GST-HRS-BF-12": 12,
                "QUE-GST-HRS-BF-15": 15,
                "QUE-GST-HRS-BF-20": 20,
                "QUE-GST-HRS-BF-25": 25,
                "QUE-GST-HRS-BF-30": 30,
                "QUE-GST-HRS-BF-35": 35,
                "QUE-GST-HRS-BF-40": 40,
                "QUE-GST-HRS-BF-16": 16,
                "QUE-GST-HRS-BF-18": 18,
                "QUE-GST-HRS-BF-22": 22,
                "QUE-GST-HRS-BF-45": 45,
                "QUE-GST-HRS-BF-F-30": 30,
                "QUE-GST-HRS-BF-F-35": 35,
                "QUE-GST-HRS-BF-F-40": 40,
            }

            body_fat = consumer_internal_markers.get("QUE-GST-HRS-BF")
            body_fat = body_fat_dict.get(body_fat)
            if gender == "M":

                if body_fat > 30:
                    pe_score += 15
                if 20 < body_fat < 30:
                    pe_score += 10
                if 12 < body_fat < 20:
                    pe_score += 5
                if 8 < body_fat < 12:
                    pe_score += 0
                if 5 < body_fat < 8:
                    pe_score += 5
                if body_fat < 5:
                    pe_score += 10
            else:
                if body_fat > 40:
                    pe_score += 15
                if 30 < body_fat < 40:
                    pe_score += 10
                if 22 < body_fat < 30:
                    pe_score += 5
                if 18 < body_fat < 22:
                    pe_score += 0
                if 15 < body_fat < 18:
                    pe_score += 5
                if body_fat < 15:
                    pe_score += 10
            if consumer_internal_markers.get("QUE-GST-HRS-EXC-Y"):
                pe_score += 5

            pe_score = (pe_score/45)/2
            total_pe = (endurence_score + pe_score)*100
            total_pe_score = 100 - total_pe
            nutrition_score.append(total_pe_score)

        if goal == ["QUE-GST-NED-PCN"]:

            pcn_score = 0

            if age > 45:
                pcn_score += 1

            severity_covid_infection = {
                "QUE-GST-RLSDSHO-RLS-SEV-SEV-M": 2,
                "QUE-GST-RLSDSHO-RLS-SEV-SEV-S": 1,
                "QUE-GST-RLSDSHO-RLS-SEV-SEV-V": 0
            }
            if consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SEV") != None:
                pcn_score += severity_covid_infection.get(
                    consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SEV"), 0)

            if bmi > 25:
                pcn_score += 1

            pcn_score = (pcn_score/4)/5

            # SEVERITY_SCALE
            severity_score = 0
            severity_scale_list_severe = ["QUE-GST-RLSDSHO-RLS-SEV-SYM-B-SE-S", "QUE-GST-RLSDSHO-RLS-SEV-SYM-F-SE-S",
                                          "QUE-GST-RLSDSHO-RLS-SEV-SYM-N-S", "QUE-GST-RLSDSHO-RLS-SEV-SYM-O-S"]
            if consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-B") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-F") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-N") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-O") in severity_scale_list_severe:
                severity_score += 3
            else:
                severity_scale_list_mild_moderate = ["QUE-GST-RLSDSHO-RLS-SEV-SYM-B-SE-MO", "QUE-GST-RLSDSHO-RLS-SEV-SYM-B-SE-M", "QUE-GST-RLSDSHO-RLS-SEV-SYM-F-SE-MO",
                                                     "QUE-GST-RLSDSHO-RLS-SEV-SYM-F-SE-M", "QUE-GST-RLSDSHO-RLS-SEV-SYM-N-MO", "QUE-GST-RLSDSHO-RLS-SEV-SYM-N-M", "QUE-GST-RLSDSHO-RLS-SEV-SYM-O-MO", "QUE-GST-RLSDSHO-RLS-SEV-SYM-O-M"]
                if consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-B") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-F") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-N") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-O") in severity_scale_list_mild_moderate:
                    severity_score += 2

            if severity_score == 0:
                severity_scale_list_no = ["QUE-GST-RLSDSHO-RLS-SEV-SYM-B-SE-N", "QUE-GST-RLSDSHO-RLS-SEV-SYM-F-SE-N",
                                          "QUE-GST-RLSDSHO-RLS-SEV-SYM-N-N", "QUE-GST-RLSDSHO-RLS-SEV-SYM-O-N"]
                if consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-B") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-F") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-N") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-O") in severity_scale_list_no:
                    severity_score += 1
            severity_score = (severity_score/3)/5

            # FREQUENCY_SCALE
            frequency_score = 0
            frequency_scale_list_vf = ["QUE-GST-RLSDSHO-RLS-SEV-SYM-C-VF", "QUE-GST-RLSDSHO-RLS-SEV-SYM-H-VF", "QUE-GST-RLSDSHO-RLS-SEV-SYM-J-VF",
                                       "QUE-GST-RLSDSHO-RLS-SEV-SYM-K-VF", "QUE-GST-RLSDSHO-RLS-SEV-SYM-L-VF", "QUE-GST-RLSDSHO-RLS-SEV-SYM-M-VF"]
            if consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-C") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-H") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-K") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-L") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-M") in frequency_scale_list_vf:
                frequency_score += 3

            else:
                frequency_scale_list_rare_freq = ["QUE-GST-RLSDSHO-RLS-SEV-SYM-C-F", "QUE-GST-RLSDSHO-RLS-SEV-SYM-C-R", "QUE-GST-RLSDSHO-RLS-SEV-SYM-H-F", "QUE-GST-RLSDSHO-RLS-SEV-SYM-H-R", "QUE-GST-RLSDSHO-RLS-SEV-SYM-J-F", "QUE-GST-RLSDSHO-RLS-SEV-SYM-J-R",
                                                  "QUE-GST-RLSDSHO-RLS-SEV-SYM-K-F", "QUE-GST-RLSDSHO-RLS-SEV-SYM-K-R", "QUE-GST-RLSDSHO-RLS-SEV-SYM-L-F", "QUE-GST-RLSDSHO-RLS-SEV-SYM-L-R", "QUE-GST-RLSDSHO-RLS-SEV-SYM-M-F", "QUE-GST-RLSDSHO-RLS-SEV-SYM-M-R"]
                if consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-C") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-H") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-K") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-L") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-M") in frequency_scale_list_rare_freq:
                    frequency_score += 2
            if frequency_score == 0:
                frequency_scale_list_no = ["QUE-GST-RLSDSHO-RLS-SEV-SYM-C-N", "QUE-GST-RLSDSHO-RLS-SEV-SYM-H-N", "QUE-GST-RLSDSHO-RLS-SEV-SYM-J-N",
                                           "QUE-GST-RLSDSHO-RLS-SEV-SYM-K-N", "QUE-GST-RLSDSHO-RLS-SEV-SYM-L-N", "QUE-GST-RLSDSHO-RLS-SEV-SYM-M-N"]
                if consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-C") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-H") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-K") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-L") or consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-M") in frequency_scale_list_no:
                    frequency_score += 1
            frequency_score = (frequency_score/3)/5

            # YES OR NO SCALE
            yes_or_no_scale = 0
            if consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-P") == "QUE-GST-RLSDSHO-RLS-SEV-SYM-P-Y":
                yes_or_no_scale += 1
            if consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-Q") == "QUE-GST-RLSDSHO-RLS-SEV-SYM-Q-Y":
                yes_or_no_scale += 1
            if consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-G") == "QUE-GST-RLSDSHO-RLS-SEV-SYM-G-Y":
                yes_or_no_scale += 1
            yes_or_no_scale = (yes_or_no_scale/3)/5

            # lung capacity
            lung_score = 0
            if consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV-SYM-R-Y"):
                lung_score += 1
            lung_score = lung_score/5
            total_pcn = (pcn_score + severity_score +
                         frequency_score + yes_or_no_scale + lung_score) * 100
            total_pcn_per = 100-total_pcn
            nutrition_score.append(total_pcn_per)

        if goal == ["QUE-GST-NED-WL"]:
            wl_score_A = 0
            current_weight = consumer_internal_markers.get("QUE-GST-CWG-KG")
            old_weight = consumer_internal_markers.get("QUE-GST-WCH-2MB-KG")
            if old_weight == None:
                old_weight = 0
            if old_weight > current_weight:
                wl_score_A += 1
            weight_change_reasons = {
                "QUE-GST-WCH-RES-DC": 4,
                "QUE-GST-WCH-RES-DK": 0,
                "QUE-GST-WCH-RES-HO": 4,
                "QUE-GST-WCH-RES-HP": 4,
                "QUE-GST-WCH-RES-MD": 4,
                "QUE-GST-WCH-RES-ME": 4,
                "QUE-GST-WCH-RES-PA": 4,
                "QUE-GST-WCH-RES-PC": 4,
                "QUE-GST-WCH-RES-PG": 4,
                "QUE-GST-WCH-RES-SD": 4,
                "QUE-GST-WCH-RES-SL": 4
            }

            if consumer_internal_markers.get("QUE-GST-WCH-RES") != None:
                wl_score_A += weight_change_reasons.get(
                    consumer_internal_markers.get("QUE-GST-WCH-RES"), 0)

            if consumer_internal_markers.get("QUE-GST-FAS-Y"):
                wl_score_A += 2

            body_type_dict = {
                "QUE-DET-BOT-ECT": 0,
                "QUE-DET-BOT-MES": 1,
                "QUE-DET-BOT-END": 2
            }
            if consumer_internal_markers.get("QUE-DET-BOT") != None:
                wl_score_A += body_type_dict.get(
                    consumer_internal_markers.get("QUE-DET-BOT"), 0)

            wl_score_A = (wl_score_A/7)/2

            height = consumer_internal_markers.get("QUE-GST-HGT-CM")/100
            ibw = 22 * height * height

            percent_ibw = ((current_weight - ibw)/ibw)/2
            weight_loss = (wl_score_A + percent_ibw)*100

            total_wl_score = 100 - weight_loss
            nutrition_score.append(round(total_wl_score, 0))

        if goal == ["QUE-GST-NED-SK"]:
            skin_score = 0
            skin_score_b = 0

            # 15 12 pts max

            if consumer_internal_markers.get("QUE-GST-SK-DL"):
                skin_score += 1
            if consumer_internal_markers.get("QUE-GST-SK-FL"):
                skin_score += 1
            if consumer_internal_markers.get("QUE-GST-SK-DLN"):
                skin_score += 1
            if consumer_internal_markers.get("QUE-GST-SK-SG"):
                skin_score += 1
            if consumer_internal_markers.get("QUE-GST-SK-AC"):
                skin_score += 1
            if consumer_internal_markers.get("QUE-GST-SK-SPO"):
                skin_score += 1
            if consumer_internal_markers.get("QUE-GST-SK-HYP"):
                skin_score += 1
            if consumer_internal_markers.get("QUE-GST-SK-HYPO"):
                skin_score += 1
            if consumer_internal_markers.get("QUE-GST-SK-PHT"):
                skin_score += 1
            if consumer_internal_markers.get("QUE-GST-SK-SNB"):
                skin_score += 1
            if consumer_internal_markers.get("QUE-GST-SK-FRC"):
                skin_score += 1
            if consumer_internal_markers.get("QUE-GST-SK-SNS"):
                skin_score += 1

            # 15.b 2 pts max
            medicine_dict = {
                "QUE-GST-MEDI-ORM": 2,
                "QUE-GST-MEDI-TPM": 2,
                "QUE-GST-MEDI-BT": 2,
                "QUE-GST-MEDI-NO": 0
            }
            if consumer_internal_markers.get("QUE-GST-MEDI") != None:
                skin_score += medicine_dict.get(
                    consumer_internal_markers.get("QUE-GST-MEDI"), 0)

            # 17 2 pts max
            sun_exposure_dict = {
                "QUE-GST-SN-EXP-LOW": 0,
                "QUE-GST-SN-EXP-MOD": 1,
                "QUE-GST-SN-EXP-HIG": 2
            }
            if consumer_internal_markers.get("QUE-GST-SN-EXP") != None:
                skin_score += sun_exposure_dict.get(
                    consumer_internal_markers.get("QUE-GST-SN-EXP"), 0)

            # 15 pts max
            if consumer_internal_markers.get("QUE-GST-SK-CON-PSO"):
                skin_score_b += 5
            if consumer_internal_markers.get("QUE-GST-SK-CON-DER"):
                skin_score_b += 5
            if consumer_internal_markers.get("QUE-GST-SK-CON-VIT"):
                skin_score_b += 5

            total_skin = (skin_score/16)*100
            total_skin = (skin_score_b + total_skin)*3
            total_skin_score = 100 - total_skin
            nutrition_score.append(total_skin_score)

        if goal == ["QUE-GST-NED-SP"] or goal == ["QUE-GST-NED-SW"]:
            sleep_score_A = 0

            # 15i max points-4
            if consumer_internal_markers.get("QUE-GST-SLS"):
                sleep_score_A += 1
            if consumer_internal_markers.get("QUE-GST-SNR"):
                sleep_score_A += 1
            if consumer_internal_markers.get("QUE-GST-SHF"):
                sleep_score_A += 1
            if consumer_internal_markers.get("QUE-GST-RLS"):
                sleep_score_A += 1

            # 15ii max score 4 points
            if consumer_internal_markers.get("QUE-GST-DYS"):
                sleep_score_A += 1
            if consumer_internal_markers.get("QUE-GST-FRZ"):
                sleep_score_A += 1
            if consumer_internal_markers.get("QUE-GST-FAPP"):
                sleep_score_A += 1
            if consumer_internal_markers.get("QUE-GST-MDF"):
                sleep_score_A += 1

            # 15iii max score 1
            if consumer_internal_markers.get("QUE-GST-HPR-Y-"):
                sleep_score_A += 1

            # 15- iv max points 3
            sleep_quality_dict = {
                "QUE-GST-SLE-VP": 3,
                "QUE-GST-SLE-PO": 2,
                "QUE-GST-SLE-FA": 1,
                "QUE-GST-SLE-GD": 0,
                "QUE-GST-SLE-VG": 0
            }
            if consumer_internal_markers.get("QUE-GST-SLE") != None:
                sleep_score_A += sleep_quality_dict.get(
                    consumer_internal_markers.get("QUE-GST-SLE"), 0)

            if caffine is not None:
                if caffine[0] > 400:  # caffine intake calculated max 4 pts
                    sleep_score_A += 4

            # step1 answer in sleep
            sleep_score_A_per = (sleep_score_A/16) * 0.25  
            insomnia_score = 0
            insomnia_dict = {
                "QUE-GST-SLPD-N": 0,
                "QUE-GST-SLPD-M": 1,
                "QUE-GST-SLPD-MO": 2,
                "QUE-GST-SLPD-S": 3,
                "QUE-GST-SLPD-VS": 4,

                "QUE-GST-SLPSA-N": 0,
                "QUE-GST-SLPSA-M": 1,
                "QUE-GST-SLPSA-MO": 2,
                "QUE-GST-SLPSA-S": 3,
                "QUE-GST-SLPSA-VS": 4,

                "QUE-GST-WKE-N": 0,
                "QUE-GST-WKE-M": 1,
                "QUE-GST-WKE-MO": 2,
                "QUE-GST-WKE-S": 3,
                "QUE-GST-WKE-VS": 4,

                "QUE-GST-SLSSAT-M": 1,
                "QUE-GST-SLSSAT-MO": 2,
                "QUE-GST-SLSSAT-N": 0,
                "QUE-GST-SLSSAT-S": 3,
                "QUE-GST-SLSSAT-VS": 4,

                "QUE-GST-SLSNOT-M": 1,
                "QUE-GST-SLSNOT-MO": 2,
                "QUE-GST-SLSNOT-N": 0,
                "QUE-GST-SLSNOT-S": 3,
                "QUE-GST-SLSNOT-VS": 4,

                "QUE-GST-SLPWOR-M": 1,
                "QUE-GST-SLPWOR-MO": 2,
                "QUE-GST-SLPWOR-N": 0,
                "QUE-GST-SLPWOR-S": 3,
                "QUE-GST-SLPWOR-VS": 4,

                "QUE-GST-SLPINT-M": 1,
                "QUE-GST-SLPINT-MO": 2,
                "QUE-GST-SLPINT-N": 0,
                "QUE-GST-SLPINT-S": 3,
                "QUE-GST-SLPINT-VS": 4,

            }

            # insomnia max points 28
            insomnia_score += insomnia_dict.get(
                consumer_internal_markers.get("QUE-GST-SLPD"), 0)
            insomnia_score += insomnia_dict.get(
                consumer_internal_markers.get("QUE-GST-SLPSA"), 0)
            insomnia_score += insomnia_dict.get(
                consumer_internal_markers.get("QUE-GST-WKE"), 0)
            insomnia_score += insomnia_dict.get(
                consumer_internal_markers.get("QUE-GST-SLSSAT"), 0)
            insomnia_score += insomnia_dict.get(
                consumer_internal_markers.get("QUE-GST-SLSNOT"), 0)
            insomnia_score += insomnia_dict.get(
                consumer_internal_markers.get("QUE-GST-SLPWOR"), 0)
            insomnia_score += insomnia_dict.get(
                consumer_internal_markers.get("QUE-GST-SLPINT"), 0)
            insomnia_score = (insomnia_score/28)/4

            # 17 apnea max points 8
            apnea_score = 0
            apnea_dict = {
                "QUE-GST-SNRLOU-Y-": 1,
                "QUE-GST-SLPTIR-Y-": 1,
                "QUE-GST-SLPSTP-Y-": 1,
                "QUE-GST-SLPBLD-Y-": 1,
                "QUE-GST-SLPNCK-Y-": 1

            }
            apnea_score += apnea_dict.get(
                consumer_internal_markers.get("QUE-GST-SNRLOU"), 0)
            apnea_score += apnea_dict.get(
                consumer_internal_markers.get("QUE-GST-SLPTIR"), 0)
            apnea_score += apnea_dict.get(
                consumer_internal_markers.get("QUE-GST-SLPSTP"), 0)
            apnea_score += apnea_dict.get(
                consumer_internal_markers.get("QUE-GST-SLPBLD"), 0)
            apnea_score += apnea_dict.get(
                consumer_internal_markers.get("QUE-GST-SLPNCK"), 0)

            if gender == "M":
                apnea_score += 1
            if bmi > 35:
                apnea_score += 1
            if consumer_internal_markers.get("QUE-GIN-AGE") > 50:
                apnea_score += 1

            apnea_score = (apnea_score/8)/4

            # rls max points 40
            rls_score = 0
            rls_dict = {

                "QUE-GST-RLSAMS-M": 2,
                "QUE-GST-RLSAMS-MD": 1,
                "QUE-GST-RLSAMS-N-": 0,
                "QUE-GST-RLSAMS-S": 3,
                "QUE-GST-RLSAMS-V": 4,

                "QUE-GST-RLSAMS-MOV-M": 1,
                "QUE-GST-RLSAMS-MOV-MD": 1,
                "QUE-GST-RLSAMS-MOV-N-": 0,
                "QUE-GST-RLSAMS-MOV-S": 3,
                "QUE-GST-RLSAMS-MOV-V": 4,

                "QUE-GST-RLSDS-MOV-M": 2,
                "QUE-GST-RLSDS-MOV-MD": 1,
                "QUE-GST-RLSDS-MOV-N-": 0,
                "QUE-GST-RLSDS-MOV-S": 3,
                "QUE-GST-RLSDS-MOV-V": 4,

                "QUE-GST-RLSDSHO-MOV-M": 2,
                "QUE-GST-RLSDSHO-MOV-MD": 1,
                "QUE-GST-RLSDSHO-MOV-N-": 0,
                "QUE-GST-RLSDSHO-MOV-S": 3,
                "QUE-GST-RLSDSHO-MOV-V": 4,

                "QUE-GST-RLSDSHO-OV-M": 2,
                "QUE-GST-RLSDSHO-OV-MD": 1,
                "QUE-GST-RLSDSHO-OV-N-": 0,
                "QUE-GST-RLSDSHO-OV-S": 3,
                "QUE-GST-RLSDSHO-OV-V": 4,


                "QUE-GST-RLSDSHO-RLS-M": 2,
                "QUE-GST-RLSDSHO-RLS-MD": 1,
                "QUE-GST-RLSDSHO-RLS-N-": 0,
                "QUE-GST-RLSDSHO-RLS-S": 3,
                "QUE-GST-RLSDSHO-RLS-V": 4,

                "QUE-GST-RLSDSHO-RLS-SEV-M": 2,
                "QUE-GST-RLSDSHO-RLS-SEV-MD": 1,
                "QUE-GST-RLSDSHO-RLS-SEV-N-": 0,
                "QUE-GST-RLSDSHO-RLS-SEV-S": 3,
                "QUE-GST-RLSDSHO-RLS-SEV-V": 4,


                "QUE-GST-RLSDSHO-RLS-SEV-IMP-M": 2,
                "QUE-GST-RLSDSHO-RLS-SEV-IMP-MD": 1,
                "QUE-GST-RLSDSHO-RLS-SEV-IMP-N-": 0,
                "QUE-GST-RLSDSHO-RLS-SEV-IMP-S": 3,
                "QUE-GST-RLSDSHO-RLS-SEV-IMP-V": 4,

                "QUE-GST-RLSDSHO-RLS-SEV-MDA-M": 2,
                "QUE-GST-RLSDSHO-RLS-SEV-MDA-MD": 1,
                "QUE-GST-RLSDSHO-RLS-SEV-MDA-N-": 0,
                "QUE-GST-RLSDSHO-RLS-SEV-MDA-S": 3,
                "QUE-GST-RLSDSHO-RLS-SEV-MDA-V": 4,

                "QUE-GST-RLSDSHO-RLS-SEV-RLS-M": 2,
                "QUE-GST-RLSDSHO-RLS-SEV-RLS-MD": 1,
                "QUE-GST-RLSDSHO-RLS-SEV-RLS-N-": 0,
                "QUE-GST-RLSDSHO-RLS-SEV-RLS-S": 3,
                "QUE-GST-RLSDSHO-RLS-SEV-RLS-V": 4,


            }
            rls_score += rls_dict.get(
                consumer_internal_markers.get("QUE-GST-RLSAMS"), 0)
            rls_score += rls_dict.get(
                consumer_internal_markers.get("QUE-GST-RLSAMS-MOV"), 0)
            rls_score += rls_dict.get(
                consumer_internal_markers.get("QUE-GST-RLSDS"), 0)
            rls_score += rls_dict.get(
                consumer_internal_markers.get("QUE-GST-RLSDSHO"), 0)
            rls_score += rls_dict.get(
                consumer_internal_markers.get("QUE-GST-RLSDSHO-OV"), 0)
            rls_score += rls_dict.get(
                consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS"), 0)
            rls_score += rls_dict.get(
                consumer_internal_markers.get("QUE-GST-RLSDSHO-RLS-SEV"), 0)
            rls_score += rls_dict.get(consumer_internal_markers.get(
                "QUE-GST-RLSDSHO-RLS-SEV-IMP"), 0)
            rls_score += rls_dict.get(consumer_internal_markers.get(
                "QUE-GST-RLSDSHO-RLS-SEV-MDA"), 0)
            rls_score += rls_dict.get(consumer_internal_markers.get(
                "QUE-GST-RLSDSHO-RLS-SEV-RLS"), 0)

            rls_score = (rls_score/40)/4

            # shiftworker max 1 pts
            shift_score = 0
            if consumer_internal_markers.get("QUE-GST-SLPSHF-Y-"):
                shift_score += 1
            shift_score = (shift_score/1)/4

            sum_score_conditions = (
                insomnia_score + apnea_score + rls_score + shift_score)*0.25  # step 2 answer
            spe_score = 0
            if spe >= 85:
                spe_score += 0
            else:
                spe_score += 5
            spe_score = (spe_score/5)/2
            sp_duration = 0
            if 7 <= sleep_duration <= 9:
                sp_duration += 0
            else:
                sp_duration += 5
            sp_duration = (sp_duration/5)/2
            sum_score_sp = (spe_score + sp_duration)*0.50  # step 3 answer
            total_sleep = (sleep_score_A_per + sum_score_sp +
                           sum_score_conditions)*100
            total_sleep_score = 100 - total_sleep
            print(total_sleep_score)
            nutrition_score.append(round(total_sleep_score))

        # gastro fitness is not yet been given so it is zero for now
        if goal == ["QUE-GST-NED-GF"]:
            gastro_fitness_score = 0
            nutrition_score.append(gastro_fitness_score)
        return nutrition_score

    @staticmethod
    def get_goal_name(goal):
        goal_dict = {
            "QUE-GST-NED-GH": "General Health",
            "QUE-GST-NED-BR": "Brain Health",
            "QUE-GST-NED-GF": "Gastro fitness",
            "QUE-GST-NED-PCN": "Post-Covid Nutrition",
            "QUE-GST-NED-PE": "Performance and Endurance",
            "QUE-GST-NED-SK": "Skin",
            "QUE-GST-NED-SP": "Sleep",
            "QUE-GST-NED-SW": "Shift Workers",
            "QUE-GST-NED-WL": "Weight management",
        }
        goal_list = [k for k, v in goal.items() if v == True]
        goal = goal_dict.get(goal_list[0])
        return goal

    # @staticmethod
    # def get_alc_quantity(harliquor_quantity,beer_quantity,wine_quantity):
    #     if harliquor_quantity ==None:
    #         harliquor_quantity = 0
    #     if beer_quantity ==None:
    #         beer_quantity = 0
    #     if wine_quantity ==None:
    #         wine_quantity = 0
    #     return ((beer_quantity/8)+(wine_quantity/8)+(harliquor_quantity/8))/30

    @staticmethod
    def get_hydration_goal(exercise_data, height):
        exercise_total_mins = sum(exercise_data.get("QUE-HEH-EAL-VO-M"))
        ml_per_cup = 240
        if exercise_total_mins:
            water_in_ml = ((height*height)*694 + \
                ((exercise_total_mins*1000)/60))
            water_in_cups = round(water_in_ml/ml_per_cup)
        else:
            water_in_ml = (height*height)*694
            water_in_cups = round(water_in_ml/ml_per_cup)

        return water_in_cups, round(water_in_ml,2)

    @staticmethod
    def cardiovascular_risk(eth, diab, gender, age, smoke, sbp, cholesterol):
        smoke = True if "QUE-GST-SMK-Y" in [k for k,
                                            v in smoke.items() if v == True] else False
        sbp = ["QUE-GST-SBP-DN"] if len([k for k, v in sbp.items() if v == True]) == 0 else [
            k for k, v in sbp.items() if v == True]
        cholesterol = ["QUE-GST-CHL-DN"]if len([k for k, v in cholesterol.items(
        ) if v == True]) == 0 else [k for k, v in cholesterol.items() if v == True]
        ethnicity = [k for k, v in eth.items() if v == True]
        cvd = 0

        sbp_const_val_dict = {
            "QUE-GST-SBP-GE80": ">=180",
            "QUE-GST-SBP-BET-160-179": "160-179",
            "QUE-GST-SBP-BET-140-159": "140-159",
            "QUE-GST-SBP-BET-120-139": "120-139",
            "QUE-GST-SBP-BET-90-120": "<120",
            "QUE-GST-SBP-LT-90": '<120',
            "QUE-GST-SBP-DN": "120-139"
        }

        sbp_selected_value = sbp_const_val_dict[sbp[0]]

        cholestrol_const_val_dict = {
            "QUE-GST-CHL-LT4": "<4",
            "QUE-GST-CHL-BET-4-4.9": "4-4.9",
            "QUE-GST-CHL-BET-5-5.9": "5-5.9",
            "QUE-GST-CHL-BET-6-6.9": "6-6.9",
            "QUE-GST-CHL-GT7": ">=7",
            "QUE-GST-CHL-DN": "4-4.9",
        }
        cholestrol_selected_value = cholestrol_const_val_dict[cholesterol[0]]
        ethnicity_const_val_dict = {
            "QUE-GIN-ETH-SEA": "south_asia",
            "QUE-GIN-ETH-EAS": "high_income_asia_pacific",
            "QUE-GIN-ETH-CAU": "central_europe",
            "QUE-GIN-ETH-HIS": "central_latin_america",
            "QUE-GIN-ETH-PAI": "central_latin_america",
            "QUE-GIN-ETH-NTA": "central_latin_america",
            "QUE-GIN-ETH-AFR": "north_africa_and_middle_east",
            "QUE-GIN-ETH-MID": "north_africa_and_middle_east"
        }
        ethnicity_selected_value = ethnicity_const_val_dict[ethnicity[0]]

        selected_gender = "women" if gender == 'F' else "men"
        selected_smoker = "smoker" if smoke else 'non_smoker'
        selected_diab = "with_diabetes" if diab else "without_diabetes"

        if 70 <= age <= 74:
            cvd = cvd_risk_helper(selected_gender, selected_smoker, selected_diab, "70-74",
                                  sbp_selected_value, cholestrol_selected_value, ethnicity_selected_value)
        elif 65 <= age <= 69:
            cvd = cvd_risk_helper(selected_gender, selected_smoker, selected_diab, "65-69",
                                  sbp_selected_value, cholestrol_selected_value, ethnicity_selected_value)
        elif 60 <= age <= 64:
            cvd = cvd_risk_helper(selected_gender, selected_smoker, selected_diab, "60-64",
                                  sbp_selected_value, cholestrol_selected_value, ethnicity_selected_value)
        elif 55 <= age <= 59:
            cvd = cvd_risk_helper(selected_gender, selected_smoker, selected_diab, "55-59",
                                  sbp_selected_value, cholestrol_selected_value, ethnicity_selected_value)
        elif 50 <= age <= 54:
            cvd = cvd_risk_helper(selected_gender, selected_smoker, selected_diab, "50-54",
                                  sbp_selected_value, cholestrol_selected_value, ethnicity_selected_value)
        elif 45 <= age <= 49:
            cvd = cvd_risk_helper(selected_gender, selected_smoker, selected_diab, "45-49",
                                  sbp_selected_value, cholestrol_selected_value, ethnicity_selected_value)
        elif 40 <= age <= 44:
            cvd = cvd_risk_helper(selected_gender, selected_smoker, selected_diab, "40-44",
                                  sbp_selected_value, cholestrol_selected_value, ethnicity_selected_value)
        elif age <= 40:
            cvd = cvd_risk_helper(selected_gender, selected_smoker, selected_diab, "40-44",
                                  sbp_selected_value, cholestrol_selected_value, ethnicity_selected_value)

        return cvd

    @staticmethod
    def get_microbiome_score(pk):
        microbiome_list = []
        try:
            microbiome_score = calc_microbiome.get_hgraph_data(
                pk)  # pass consumer id instead of this string
        except:
            return None
        goals = ["Physiology", "Mental Health", "Aesthetics", "Overall health", "Sleep", "Brain Health", "Skin", "Performance & Endurance", "Weight loss",
                 "Post COVID care", "Shift workers", "General health", "Gastro fitness"]
        h_graph_focus_area = microbiome_score["h_graph_focus_area"]
        for i in goals:
            for j in range(len(h_graph_focus_area)):
                if h_graph_focus_area[j]["Class"] == i:
                    microbiome_list.append(
                        h_graph_focus_area[j]["Overall_health_score"])
        return microbiome_list

    @staticmethod
    def get_genetic_score(pk):
        genetic_score_list = []
        goals = ["Physiology", "Mental_Health", "Aesthetics", "Overall_Health", "Sleep", "Brain Health", "Skin", "Performance & Endurance", "Weight loss",
                 "Post COVID care", "Shift workers", "General health", "Gastro fitness"]
        try:
            genetic_score = calc_genetics.get_overall_health_data(
                pk)  # Need to pass the consumer_id here
            focus_area = genetic_score[0]
            for i in goals:
                if i in focus_area.keys():
                    genetic_score_list.append(focus_area[i])
            genetic_score_list.append(focus_area["Goal_Score"])
            return genetic_score_list
        except:
            return None
        

    @staticmethod
    def get_microbiome_disease_risk(pk):  # used for suppliment recommentation
        li = ["Vitamin B1", "Vitamin B2",
              "Vitamin B5", "Vitamin B6", "Vitamin B9"]
        microbiome_dict = {}
        try:
            microbiome_score = calc_microbiome.get_disease_risk_data(
                pk)
            h_graph_focus_area = microbiome_score["h_graph_focus_area"]
            for i in li:
                for j in range(len(h_graph_focus_area)):
                    if h_graph_focus_area[j]["Phenotype"] == i:
                        microbiome_dict[i] = h_graph_focus_area[j]["Result"]
            return microbiome_dict
              # pass consumer id instead of this string
        except:
            return microbiome_dict
        

    @staticmethod
    def get_genetics_disease_risk(pk):  # used for suppliment recommentation
        li = ["Vitamin B9", "Vitamin B6", "Vitamin B12", "Vitamin A",
              "Vitamin C", "Vitamin D", "Vitamin E", "Vitamin B2", "Iron"]
        genetic_dict = {}
        try:
            genetic_score = calc_genetics.get_genetics_risk_area(
                pk)  # pass consumer id instead of this string
            l = genetic_score["overall_health_score"]
            for i in li:
                for j in range(len(l)):
                    if l[j]["Condition"] == i:
                        if l[j]["Color"] == "Red":
                            genetic_dict[i] = "low"
                        elif l[j]["Color"] == "Orange":
                            genetic_dict[i] = "moderate"
                        elif l[j]["Color"] == "Green":
                            genetic_dict[i] = "high"
            return genetic_dict
        except:
            return genetic_dict

        

    @staticmethod
    def get_diet_type(goal_data, age, gender):
        goal = [k for k, v in goal_data.items() if v == True]
        diet_type = None

        if goal == ["QUE-GST-NED-GH"] and age > 18 and gender == "M":
            diet_type = diet_rules.general_health[0]
        else:
            if goal == ["QUE-GST-NED-GH"] and gender == "F":
                diet_type = diet_rules.general_health[1]

        if goal == ["QUE-GST-NED-BR"]:  # need work
            diet_type = diet_rules.brain_health[0]

        elif goal == ["QUE-GST-NED-GF"]:
            diet_type = diet_rules.gastro_fitness[0]

        elif goal == ["QUE-GST-NED-PCN"]:
            diet_type = diet_rules.post_covid_nutrition[0]

        elif goal == ["QUE-GST-NED-PE"]:
            diet_type = diet_rules.performance_and_endurance[0]

        elif goal == ["QUE-GST-NED-SK"]:
            diet_type = diet_rules.skin[0]

        elif goal == ["QUE-GST-NED-SP"]:
            diet_type = diet_rules.sleep[0]  # needs work

        elif goal == ["QUE-GST-NED-SW"]:
            diet_type = diet_rules.general_health[0]

        elif goal == ["QUE-GST-NED-WL"] and age > 18:
            diet_type = diet_rules.weight_loss[1]
        else:
            if goal == ["QUE-GST-NED-WL"]:
                diet_type = diet_rules.weight_loss[0]

        return diet_type

    @staticmethod
    def suppliment_recommendations(age, pregnant, lactating, smoker, rda_met_percent, gender, category, genetics=None, microbiome=None):
        rda_male = {
            "Vit A": 900.0,
            "Vit B1": 1.2,
            "Vit B2": 1.3,
            "Vit B3": 16.0,
            "Vit B5": 5.0,
            "Vit B6": 1.3,
            "Vit B7-Biotin": 30.0,
            "Vit B9-Folic acid": 400.0,
            "Vit B12": 2.4,
            "Vit C": 125.0,
            "Vit D": 15.0,
            "Vit E": 15.0,
            "Vit K": 120.0,
            "Iron": 8.0,
            "Zinc": 11.0,
            "Magnesium": 420.0,
            "Calcium": 1000.0,
            "Selenium": 55.0,
            "Copper": 900.0,
            "Manganese": 2.3,
            "Chromium": 35.0,
            "Phosphorus": 700.0,
        }
        rda_female = {
            "Vit A": 700.0,
            "Vit B1": 1.1,
            "Vit B2": 1.3,
            "Vit B3": 14.0,
            "Vit B5": 5.0,
            "Vit B6": 1.3,
            "Vit B7-Biotin": 30.0,
            "Vit B9-Folic acid": 400.0,
            "Vit B12": 2.4,
            "Vit C": 75.0,
            "Vit D": 15.0,
            "Vit E": 15.0,
            "Vit K": 90.0,
            "Iron": 18.0,
            "Zinc": 11.0,
            "Magnesium": 320.0,
            "Calcium": 1000.0,
            "Selenium": 55.0,
            "Copper": 900.0,
            "Manganese": 1.8,
            "Chromium": 25.0,
            "Phosphorus": 700.0,
        }

        gender = "male" if gender == "M" else "female"
        if gender == "male":
            rda = rda_male.get(category, {})
        elif gender == "female":
            rda = rda_female.get(category, {})
        # rda = Calculate.get_rda(
        #     age, rda_constants, pregnant, lactating, smoker)
        # rda = constants.rda.get(gender).get(category, 0)
        unit = constants.suppliment_unit.get(category, 0)
        tul = constants.tul.get(category, 0)
        bio_availability = constants.bio_percent.get(category, 0)
        formula = 0
        risk_genetics = {
            "low": 0,
            "high": 10,
            "moderate": 5,
            None: 0
        }
        risk_microbiome = {
            "low-risk": 0,
            "moderate-risk": 5,
            "high-risk": 10,
            None: 0
        }
        replenishment_dose = ((rda*(100-rda_met_percent)-(
            risk_genetics[genetics] + risk_microbiome[microbiome]))/100)*(100/bio_availability)
        maintanence_dose = ((rda*(100-rda_met_percent)-(
            risk_genetics[genetics] + risk_microbiome[microbiome]))/100)*(100/bio_availability)
        if 0 <= rda_met_percent <= 30:
            formula = tul

        elif 30 < rda_met_percent <= 50:
            formula = replenishment_dose + maintanence_dose
        elif 50 < rda_met_percent <= 75:
            formula = replenishment_dose + (maintanence_dose/2)
        elif 75 < rda_met_percent < 100:
            formula = replenishment_dose + (maintanence_dose/2)

        elif rda_met_percent >= 100:
            formula = replenishment_dose + (maintanence_dose/2)

        if formula < 0:
            formula = 0

        if formula > tul:
            if category == "Vit C":
                formula = tul/2
            else:
                formula = tul

        if category == "Calcium":
            formula = formula/2
        formula=format(formula, '.2f')    
        return formula, unit, category

    @staticmethod
    def get_suppliment_recommendation_pack_2(goal):

        strain_1 = constants.strain_1.get(goal)
        dose = constants.dose.get(goal)

        # nutrition scoring to get secondary goal

        return strain_1, dose

    @staticmethod
    def get_secondary_recommendation_pack_2(goal_data, consumer_internal_markers, smoke_data, grains, fat, alcohol, gender, protein, veg_and_fru_servings, pal, spe, sleep_duration, bmi, ethnicity, pk):
        # pk = "1.002.0004"
        goal =""
        try:
            microbiome_data = calc_microbiome.get_disease_risk_data(
                pk)  # pass consumer id instead of this string
            data = microbiome_data["h_graph_focus_area"]
            general_health_list = ["Gluten intolerance",
                                   "NAFLD",
                                   "Lactose intolerance",
                                   "Vitamin B9",
                                   "Vitamin B6",
                                   "Vitamin A",
                                   "Vitamin D",
                                   "Calcium",
                                   "Alcoholic liver disease",
                                   "Anorexia nervosa",
                                   "Blood Lipid levels-HDL",
                                   "Blood Lipid levels-TGL",
                                   "CVD",
                                   "Gall stone",
                                   "Graves disease",
                                   "Hypertension",
                                   "Insulin sensitivity",
                                   "Kidney Stones",
                                   "Metabolic syndrome",
                                   "Osteoporosis",
                                   "Prediabetes",
                                   "RA",
                                   "Spondyloarthritis",
                                   "Type II Diabetes",
                                   "Vitamin B1",
                                   "Vitamin B2",
                                   "Vitamin B5"
                                   ]
            skin_health_list = ["Acne",
                                "Atopic dermatitis",
                                "Psoriasis",
                                "Ageing",
                                "Vitiligo",
                                ]
            shift_workers_list = ["Sleep quality"]
            post_covid_list = ["Covid severity"]
            brain_health_list = ["Anxiety",
                                 "Attention-deficit-hyperactive disorder",
                                 "Migraine",
                                 "Depression",
                                 "PTSD",
                                 ]
            performance_list = ["Aerobic performance",
                                "Muscle strength",
                                "Active lifestyle",
                                "Endurance",
                                "Fibromyalgia",
                                ]
            weight_managenment_list = ["Obesity",
                                       "Carbohydrate metabolism",
                                       "Fat metabolism",
                                       "Protein metabolism",
                                       ]
            sleep_health_list = ["Sleep quality"]
            gastro_gerd_list = ["Irritable Bowel Syndrome",
                                "Inflammatory Bowel Disease",
                                "Bloating",
                                "Constipation",
                                "Crohns Disease",
                                "Diarrhea",
                                "Flatulence",
                                "Ulcerative Colitis",
                                "GERD",
                                ]
            risk_type = ["Poor","moderate","Good"]
            gen_health = secondary_goal_microbiome_helper(general_health_list,data,risk_type)
            skin_health = secondary_goal_microbiome_helper(skin_health_list,data,risk_type)
            shift_workers = secondary_goal_microbiome_helper(shift_workers_list,data,risk_type)
            post_covid = secondary_goal_microbiome_helper(post_covid_list,data,risk_type)
            brain_health = secondary_goal_microbiome_helper(brain_health_list,data,risk_type)
            performance = secondary_goal_microbiome_helper(performance_list,data,risk_type)
            weight_managenment = secondary_goal_microbiome_helper(weight_managenment_list,data,risk_type)
            sleep_health = secondary_goal_microbiome_helper(sleep_health_list,data,risk_type)
            gastro_gerd = secondary_goal_microbiome_helper(gastro_gerd_list,data,risk_type)

            if "high_risk" in gen_health or "high_risk" in skin_health or "high_risk" in shift_workers or "high_risk" in post_covid or "high_risk" in brain_health or "high_risk" in performance or "high_risk" in weight_managenment or "high_risk" in sleep_health or "high_risk" in gastro_gerd:
    
                if goal_data == "General Health":
                    general_health_priority_list = [
                        "Gastro fitness", "Weight management", "Performance & endurance", "Sleep","Brain Health","Skin","Shift workers"]
                    high_risk_list_count_gen_health = [gastro_gerd.count("high_risk"), weight_managenment.count(
                        "high_risk"), performance.count("high_risk"), sleep_health.count("high_risk"),brain_health.count("high_risk"),skin_health.count("high_risk"),shift_workers.count("high_risk")]
                    goal = secondary_goal_priority_helper(
                        general_health_priority_list, high_risk_list_count_gen_health)

                elif goal_data == "Weight management":
                    weight_mgmt_priority_list = [
                        "General Health", "Gastro fitness", "Performance & endurance", "Sleep","Brain Health","Skin", "Shift workers"]
                    high_risk_list_count_weight = [gen_health.count("high_risk"), gastro_gerd.count(
                        "high_risk"), performance.count("high_risk"), sleep_health.count("high_risk"),brain_health.count("high_risk"),skin_health.count("high_risk"), shift_workers.count("high_risk")]
                    goal = secondary_goal_priority_helper(
                        weight_mgmt_priority_list, high_risk_list_count_weight)
                        # ----------------------------------------------------------------------------------------------------------------------

                elif goal_data == "Shift workers":
                    shift_worker_priority_list = [
                        "Sleep","Brain Health", "General Health", "Gastro fitness", "Weight management", "Performance & endurance","Skin" ]
                    high_risk_list_count_shift = [sleep_health.count("high_risk"),brain_health.count("high_risk"), gen_health.count(
                        "high_risk"), gastro_gerd.count("high_risk"), weight_managenment.count("high_risk"), performance.count("high_risk"),skin_health.count("high_risk")]
                    goal = secondary_goal_priority_helper(
                        shift_worker_priority_list, high_risk_list_count_shift)
                    # ----------------------------------------------------------------------------------------------------------------------
                elif goal_data == "Sleep":
                    sleep_priority_list = ["Brain Health","Weight management", "Gastro fitness",
                                           "General Health","Skin", "Performance & endurance", "Shift workers"]
                    high_risk_list_count_sleep = [brain_health.count("high_risk"),weight_managenment.count("high_risk"), gastro_gerd.count(
                        "high_risk"), gen_health.count("high_risk"),skin_health.count("high_risk"), performance.count("high_risk"), shift_workers.count("high_risk")]
                    goal = secondary_goal_priority_helper(
                        sleep_priority_list, high_risk_list_count_sleep)
                    # ----------------------------------------------------------------------------------------------------------------------
                elif goal_data == "Performance & endurance":
                    performance_priority_list = [
                        "General Health", "Sleep", "Gastro fitness","Brain Health", "Weight management","Skin", "Shift workers"]
                    high_risk_list_count_performance = [gen_health.count("high_risk"), sleep_health.count(
                        "high_risk"), gastro_gerd.count("high_risk"),brain_health.count("high_risk"), weight_managenment.count("high_risk"),skin_health.count("high_risk"), shift_workers.count("high_risk")]
                    goal = secondary_goal_priority_helper(
                        performance_priority_list, high_risk_list_count_performance)
                    # ----------------------------------------------------------------------------------------------------------------------
                elif goal_data == "Gastro fitness":
                    gastro_priority_list = ["Weight management", "Sleep", "Brain Health","General Health", "Performance & endurance","Skin" "Shift workers"]
                        
                    high_risk_list_count_gastro = [weight_managenment.count("high_risk"), sleep_health.count(
                        "high_risk"), brain_health.count("high_risk"), gen_health.count("high_risk"),performance.count("high_risk"), skin_health.count("high_risk"),shift_workers.count("high_risk")]
                    goal = secondary_goal_priority_helper(
                        gastro_priority_list, high_risk_list_count_gastro)
                    # ----------------------------------------------------------------------------------------------------------------------
                elif goal_data == "Brain Health":
                    brain_priority_list = ["Sleep","General Health","Gastro fitness", "Performance & endurance", "Weight management", "Skin", "Shift workers"]
                    high_risk_list_count_brain = [ sleep_health.count("high_risk"),gen_health.count("high_risk"), gastro_gerd.count("high_risk"), performance.count(
                        "high_risk"), weight_managenment.count("high_risk"),skin_health.count("high_risk"), shift_workers.count("high_risk")]
                    goal = secondary_goal_priority_helper(
                        brain_priority_list, high_risk_list_count_brain)
                    # ----------------------------------------------------------------------------------------------------------------------
        
                elif goal_data == "Skin":
                    brain_priority_list = ["General Health","Gastro fitness","Sleep", "Weight management","Performance & endurance", "Brain Health", "Shift workers"]
                    high_risk_list_count_brain = [ gen_health.count("high_risk"), gastro_gerd.count("high_risk"), sleep_health.count("high_risk"), weight_managenment.count("high_risk"),performance.count(
                        "high_risk"),brain_health.count("high_risk"), shift_workers.count("high_risk")]
                    goal = secondary_goal_priority_helper(
                        brain_priority_list, high_risk_list_count_brain)
            # moderate risk---------------
            elif "moderate_risk" in gen_health or "moderate_risk" in skin_health or "moderate_risk" in shift_workers or "moderate_risk" in post_covid or "moderate_risk" in brain_health or "moderate_risk" in performance or "moderate_risk" in weight_managenment or "moderate_risk" in sleep_health or "moderate_risk" in gastro_gerd:
        
                if goal_data == "General Health":
                    general_health_priority_list = [
                        "Gastro fitness", "Weight management", "Performance & endurance", "Sleep","Brain Health","Skin","Shift workers"]
                    moderate_risk_list_count_gen_health = [gastro_gerd.count("moderate_risk"), weight_managenment.count(
                        "moderate_risk"), performance.count("moderate_risk"), sleep_health.count("moderate_risk"),brain_health.count("moderate_risk"),skin_health.count("moderate_risk"),shift_workers.count("moderate_risk")]
                    goal = secondary_goal_priority_helper(
                        general_health_priority_list, moderate_risk_list_count_gen_health)

                elif goal_data == "Weight management":
                    weight_mgmt_priority_list = [
                        "General Health", "Gastro fitness", "Performance & endurance", "Sleep","Brain Health","Skin", "Shift workers"]
                    moderate_risk_list_count_weight = [gen_health.count("moderate_risk"), gastro_gerd.count(
                        "moderate_risk"), performance.count("moderate_risk"), sleep_health.count("moderate_risk"),brain_health.count("moderate_risk"),skin_health.count("moderate_risk"), shift_workers.count("moderate_risk")]
                    goal = secondary_goal_priority_helper(
                        weight_mgmt_priority_list, moderate_risk_list_count_weight)
                        # ----------------------------------------------------------------------------------------------------------------------

                elif goal_data == "Shift workers":
                    shift_worker_priority_list = [
                        "Sleep","Brain Health", "General Health", "Gastro fitness", "Weight management", "Performance & endurance","Skin" ]
                    moderate_risk_list_count_shift = [sleep_health.count("moderate_risk"),brain_health.count("moderate_risk"), gen_health.count(
                        "moderate_risk"), gastro_gerd.count("moderate_risk"), weight_managenment.count("moderate_risk"), performance.count("moderate_risk"),skin_health.count("moderate_risk")]
                    goal = secondary_goal_priority_helper(
                        shift_worker_priority_list, moderate_risk_list_count_shift)
                    # ----------------------------------------------------------------------------------------------------------------------
                elif goal_data == "Sleep":
                    sleep_priority_list = ["Brain Health","Weight management", "Gastro fitness",
                                           "General Health","Skin", "Performance & endurance", "Shift workers"]
                    moderate_risk_list_count_sleep = [brain_health.count("moderate_risk"),weight_managenment.count("moderate_risk"), gastro_gerd.count(
                        "moderate_risk"), gen_health.count("moderate_risk"),skin_health.count("moderate_risk"), performance.count("moderate_risk"), shift_workers.count("moderate_risk")]
                    goal = secondary_goal_priority_helper(
                        sleep_priority_list, moderate_risk_list_count_sleep)
                    # ----------------------------------------------------------------------------------------------------------------------
                elif goal_data == "Performance & endurance":
                    performance_priority_list = [
                        "General Health", "Sleep", "Gastro fitness","Brain Health", "Weight management","Skin", "Shift workers"]
                    moderate_risk_list_count_performance = [gen_health.count("moderate_risk"), sleep_health.count(
                        "moderate_risk"), gastro_gerd.count("moderate_risk"),brain_health.count("moderate_risk"), weight_managenment.count("moderate_risk"),skin_health.count("moderate_risk"), shift_workers.count("moderate_risk")]
                    goal = secondary_goal_priority_helper(
                        performance_priority_list, moderate_risk_list_count_performance)
                    # ----------------------------------------------------------------------------------------------------------------------
                elif goal_data == "Gastro fitness":
                    gastro_priority_list = ["Weight management", "Sleep", "Brain Health","General Health", "Performance & endurance","Skin" "Shift workers"]
                        
                    moderate_risk_list_count_gastro = [weight_managenment.count("moderate_risk"), sleep_health.count(
                        "moderate_risk"), brain_health.count("moderate_risk"), gen_health.count("moderate_risk"),performance.count("moderate_risk"), skin_health.count("moderate_risk"),shift_workers.count("moderate_risk")]
                    goal = secondary_goal_priority_helper(
                        gastro_priority_list, moderate_risk_list_count_gastro)
                    # ----------------------------------------------------------------------------------------------------------------------
                elif goal_data == "Brain Health":
                    brain_priority_list = ["Sleep","General Health","Gastro fitness", "Performance & endurance", "Weight management", "Skin", "Shift workers"]
                    moderate_risk_list_count_brain = [ sleep_health.count("moderate_risk"),gen_health.count("moderate_risk"), gastro_gerd.count("moderate_risk"), performance.count(
                        "moderate_risk"), weight_managenment.count("moderate_risk"),skin_health.count("moderate_risk"), shift_workers.count("moderate_risk")]
                    goal = secondary_goal_priority_helper(
                        brain_priority_list, moderate_risk_list_count_brain)
                    # ----------------------------------------------------------------------------------------------------------------------
        
                elif goal_data == "Skin":
                    brain_priority_list = ["General Health","Gastro fitness","Sleep", "Weight management","Performance & endurance", "Brain Health", "Shift workers"]
                    moderate_risk_list_count_brain = [ gen_health.count("moderate_risk"), gastro_gerd.count("moderate_risk"), sleep_health.count("moderate_risk"), weight_managenment.count("moderate_risk"),performance.count(
                        "moderate_risk"),brain_health.count("moderate_risk"), shift_workers.count("moderate_risk")]
                    goal = secondary_goal_priority_helper(
                        brain_priority_list, moderate_risk_list_count_brain)
            # low risk------------------
            elif "low_risk" in gen_health or "low_risk" in skin_health or "low_risk" in shift_workers or "low_risk" in post_covid or "low_risk" in brain_health or "low_risk" in performance or "low_risk" in weight_managenment or "low_risk" in sleep_health or "low_risk" in gastro_gerd:
            
                if goal_data == "General Health":
                    general_health_priority_list = [
                        "Gastro fitness", "Weight management", "Performance & endurance", "Sleep","Brain Health","Skin","Shift workers"]
                    low_risk_list_count_gen_health = [gastro_gerd.count("low_risk"), weight_managenment.count(
                        "low_risk"), performance.count("low_risk"), sleep_health.count("low_risk"),brain_health.count("low_risk"),skin_health.count("low_risk"),shift_workers.count("low_risk")]
                    goal = secondary_goal_priority_helper(
                        general_health_priority_list, low_risk_list_count_gen_health)

                elif goal_data == "Weight management":
                    weight_mgmt_priority_list = [
                        "General Health", "Gastro fitness", "Performance & endurance", "Sleep","Brain Health","Skin", "Shift workers"]
                    low_risk_list_count_weight = [gen_health.count("low_risk"), gastro_gerd.count(
                        "low_risk"), performance.count("low_risk"), sleep_health.count("low_risk"),brain_health.count("low_risk"),skin_health.count("low_risk"), shift_workers.count("low_risk")]
                    goal = secondary_goal_priority_helper(
                        weight_mgmt_priority_list, low_risk_list_count_weight)
                        # ----------------------------------------------------------------------------------------------------------------------

                elif goal_data == "Shift workers":
                    shift_worker_priority_list = [
                        "Sleep","Brain Health", "General Health", "Gastro fitness", "Weight management", "Performance & endurance","Skin" ]
                    low_risk_list_count_shift = [sleep_health.count("low_risk"),brain_health.count("low_risk"), gen_health.count(
                        "low_risk"), gastro_gerd.count("low_risk"), weight_managenment.count("low_risk"), performance.count("low_risk"),skin_health.count("low_risk")]
                    goal = secondary_goal_priority_helper(
                        shift_worker_priority_list, low_risk_list_count_shift)
                    # ----------------------------------------------------------------------------------------------------------------------
                elif goal_data == "Sleep":
                    sleep_priority_list = ["Brain Health","Weight management", "Gastro fitness",
                                           "General Health","Skin", "Performance & endurance", "Shift workers"]
                    low_risk_list_count_sleep = [brain_health.count("low_risk"),weight_managenment.count("low_risk"), gastro_gerd.count(
                        "low_risk"), gen_health.count("low_risk"),skin_health.count("low_risk"), performance.count("low_risk"), shift_workers.count("low_risk")]
                    goal = secondary_goal_priority_helper(
                        sleep_priority_list, low_risk_list_count_sleep)
                    # ----------------------------------------------------------------------------------------------------------------------
                elif goal_data == "Performance & endurance":
                    performance_priority_list = [
                        "General Health", "Sleep", "Gastro fitness","Brain Health", "Weight management","Skin", "Shift workers"]
                    low_risk_list_count_performance = [gen_health.count("low_risk"), sleep_health.count(
                        "low_risk"), gastro_gerd.count("low_risk"),brain_health.count("low_risk"), weight_managenment.count("low_risk"),skin_health.count("low_risk"), shift_workers.count("low_risk")]
                    goal = secondary_goal_priority_helper(
                        performance_priority_list, low_risk_list_count_performance)
                    # ----------------------------------------------------------------------------------------------------------------------
                elif goal_data == "Gastro fitness":
                    gastro_priority_list = ["Weight management", "Sleep", "Brain Health","General Health", "Performance & endurance","Skin" "Shift workers"]
                        
                    low_risk_list_count_gastro = [weight_managenment.count("low_risk"), sleep_health.count(
                        "low_risk"), brain_health.count("low_risk"), gen_health.count("low_risk"),performance.count("low_risk"), skin_health.count("low_risk"),shift_workers.count("low_risk")]
                    goal = secondary_goal_priority_helper(
                        gastro_priority_list, low_risk_list_count_gastro)
                    # ----------------------------------------------------------------------------------------------------------------------
                elif goal_data == "Brain Health":
                    brain_priority_list = ["Sleep","General Health","Gastro fitness", "Performance & endurance", "Weight management", "Skin", "Shift workers"]
                    low_risk_list_count_brain = [ sleep_health.count("low_risk"),gen_health.count("low_risk"), gastro_gerd.count("low_risk"), performance.count(
                        "low_risk"), weight_managenment.count("low_risk"),skin_health.count("low_risk"), shift_workers.count("low_risk")]
                    goal = secondary_goal_priority_helper(
                        brain_priority_list, low_risk_list_count_brain)
                    # ----------------------------------------------------------------------------------------------------------------------
        
                elif goal_data == "Skin":
                    brain_priority_list = ["General Health","Gastro fitness","Sleep", "Weight management","Performance & endurance", "Brain Health", "Shift workers"]
                    low_risk_list_count_brain = [ gen_health.count("low_risk"), gastro_gerd.count("low_risk"), sleep_health.count("low_risk"), weight_managenment.count("low_risk"),performance.count(
                        "low_risk"),brain_health.count("low_risk"), shift_workers.count("low_risk")]
                    goal = secondary_goal_priority_helper(
                        brain_priority_list, low_risk_list_count_brain)
            return goal
        
        except:
            ethnicity = [k for k, v in ethnicity.items() if v == True]
            general_health_list = []
            pe_list = []
            weight_mgmt_list = []
            sleep_list = []
            gastro_list = []
            shift_worker_list = []
            water_score=0
            smoke_ = [k for k, v in smoke_data.items() if v == True]
            smoke = {
                'QUE-GST-SMK-Y': 0,
                'QUE-GST-SMK-X': 1,
                'QUE-GST-SMK-N': 3,
                None: ""
            }
            try:
                smoke_score = smoke.get(smoke_[0])
            except:
                smoke_score = ""
            if smoke_score == 3:
                general_health_list.append("low_risk")
            elif smoke_score == 1:
                general_health_list.append("moderate_risk")
            elif smoke_score == 0:
                general_health_list.append("high_risk")
            else:
                pass
            # for general health
            grains_dict = {
                0: "high_risk",
                5: "low_risk"
            }
            grains_score = 0
            grain_list = ["QUE-HEH-MEA-BF-FP-GR","QUE-HEH-MEA-LN-FP-GR","QUE-HEH-MEA-DN-FP-GR"]
            grain_percentage = 0
            for item in grain_list:
                if consumer_internal_markers.get(item) != None:
                    grain_percentage += consumer_internal_markers.get(item)
            grains = grain_percentage / len(grain_list)
            if grains > 20.2:
                grains_score = 5
            elif grains < 20:
                grains_score = 0

            general_health_list.append(grains_dict.get(grains_score))

            # fat intake percentage max score 10 pts
            fat_dict = {
                0: "high_risk",
                10: "low_risk",
                5: "moderate_risk"
            }
            if 25 < fat < 35:
                fat_score = 10
            elif fat < 25:
                fat_score = 5
            elif fat > 35:
                fat_score = 0
            general_health_list.append(fat_dict.get(fat_score))

            # overll health alcohol consumption (daily intake/8 = drinks)
            alc_dict = {
                0: "high_risk",
                10: "low_risk",
                5: "moderate_risk"
            }
            if alcohol != None:
                alcohol = alcohol/8
                # alc_score = 0
                print(gender)
                if alcohol == 0:
                    alc_score = 10
                elif gender == "M":
                    if 1 <= alcohol <= 2:
                        alc_score = 5
                    elif alcohol > 2:
                        alc_score = 0
                    else:
                        alc_score = 0
                elif gender == "F":
                    if 0 < alcohol <= 1:
                        alc_score = 5
                    elif alcohol > 1:
                        alc_score = 0
                    else:
                        alc_score = 0
                general_health_list.append(alc_dict.get(alc_score))

            # protein overallhealth
            protein_score = 0
            protein_dict = {
                0: "high_risk",
                10: "low_risk",
                5: "moderate_risk"
            }
            protein_list = ["QUE-HEH-MEA-BF-FP-PR","QUE-HEH-MEA-LN-FP-PR","QUE-HEH-MEA-DN-FP-PR"]
            protein_percentage = 0
            for item in protein_list:
                if consumer_internal_markers.get(item) != None:
                    protein_percentage += consumer_internal_markers.get(item)
            protein = protein_percentage / len(protein_list)
            if 10 < protein < 35:
                protein_score = 10
            elif 35 > protein < 40:
                protein_score = 5
            elif protein > 40:
                protein_score = 0
            else:
                protein_score = 0
            general_health_list.append(protein_dict.get(protein_score))

            # water_consumption 0verallhealth
            water_intake = consumer_internal_markers.get(
                "QUE-HEH-LIQ-WA-QN", 0)
            water_dict = {
                0: "high_risk",
                5: "low_risk",
                3: "moderate_risk"
            }
            if water_intake != None:
                water_intake = int(water_intake)
                if gender == "F":
                    if water_intake != 0:
                        if water_intake > 2700:
                            water_score = 5
                        elif 2000 <= water_intake < 2700:
                            water_score = 3
                        elif water_intake < 2000:
                            water_score = 0
                elif gender == "M":
                    if water_intake != 0:
                        if water_intake >= 3700:
                            water_score = 5
                        elif 2500 <= water_intake < 3700:
                            water_score = 3
                        elif water_intake < 2500:
                            water_score = 0
                general_health_list.append(water_dict.get(water_score))

            # (vegetable and fruit intake/100 = servings)
            veg_dict = {
                0: "high_risk",
                10: "low_risk",
                5: "moderate_risk"
            }
            if veg_and_fru_servings > 5:
                veg_and_fru_servings_score = 10
            elif 2 <= veg_and_fru_servings < 5:
                veg_and_fru_servings_score = 5
            elif veg_and_fru_servings < 2:
                veg_and_fru_servings_score = 0
            general_health_list.append(
                veg_dict.get(veg_and_fru_servings_score))

            # weight management
            bmi_dict = {
                0: "high_risk",
                10: "low_risk",
                5: "moderate_risk"
            }
            if ethnicity == ["QUE-GIN-ETH-AFR"]:  # bmi classification
                if bmi < 21.9:
                    bmi_score = 5
                elif 22 <= bmi < 28:
                    bmi_score = 10
                elif 28 <= bmi < 33:
                    bmi_score = 5
                elif bmi > 33:
                    bmi_score = 0
                weight_mgmt_list.append(bmi_dict.get(bmi_score))

            elif ethnicity == ["QUE-GIN-ETH-CAU"]:
                if bmi < 18:
                    bmi_score = 5
                elif 18 <= bmi < 25:
                    bmi_score = 10
                elif 25 <= bmi < 30:
                    bmi_score = 5
                elif bmi > 30:
                    bmi_score = 0
                weight_mgmt_list.append(bmi_dict.get(bmi_score))

            elif ethnicity == ["QUE-GIN-ETH-EAS"]:
                if bmi < 18.5:
                    bmi_score = 5
                elif 18.5 <= bmi <= 23:
                    bmi_score = 10
                elif 23 < bmi <= 25:
                    bmi_score = 5
                elif bmi >= 25:
                    bmi_score = 0
                weight_mgmt_list.append(bmi_dict.get(bmi_score))

            elif ethnicity == ["QUE-GIN-ETH-HIS"] or ethnicity == ["QUE-GIN-ETH-NTA"] or ethnicity == ["QUE-GIN-ETH-PAI"]:
                if bmi < 18.4:
                    bmi_score = 5
                elif 18.5 <= bmi < 25:
                    bmi_score = 10
                elif 25 <= bmi < 30:
                    bmi_score = 5
                elif bmi >= 30:
                    bmi_score = 0
                weight_mgmt_list.append(bmi_dict.get(bmi_score))

            elif ethnicity == ["QUE-GIN-ETH-MID"]:
                if bmi < 18:
                    bmi_score = 5
                elif 18 <= bmi < 23:
                    bmi_score = 10
                elif 23 <= bmi <= 27.5:
                    bmi_score = 5
                elif bmi >= 27.6:
                    bmi_score = 0
                weight_mgmt_list.append(bmi_dict.get(bmi_score))

            elif ethnicity == ["QUE-GIN-ETH-SEA"]:
                if bmi < 18:
                    bmi_score = 5
                elif 18 <= bmi < 25:
                    bmi_score = 10
                elif 25 <= bmi < 30:
                    bmi_score = 5
                elif bmi >= 30:
                    bmi_score = 0
                weight_mgmt_list.append(bmi_dict.get(bmi_score))

            # for performance and endurence
            pal_dict = {
                10: "low_risk",
                5: "moderate_risk",
                0: "high_risk"
            }
            if pal < 1.20:
                pal_score = 0
            elif 1.20 <= pal < 1.375:
                pal_score = 5
            elif 1.375 <= pal < 1.55:
                pal_score = 10
            elif 1.55 <= pal < 1.725:
                pal_score = 15
            elif pal >= 1.725:
                pal_score = 20
            pe_list.append(pal_dict.get(pal_score))

            # sleep
            sleep_score_dict = {
                5: "low_risk",
                3: "moderate_risk",
                0: "high_risk"
            }
            if spe >= 85:
                sleep_score = 5
            elif 70 <= spe <= 84.99:
                sleep_score = 3
            elif spe < 70:
                sleep_score = 0
            sleep_list.append(sleep_score_dict.get(sleep_score))
            shift_worker_list.append(sleep_score_dict.get(sleep_score))

            # sleep_duration max pts 5 mental health
            duration_score_dict = {
                5: "low_risk",
                3: "moderate_risk",
                0: "high_risk"

            }
            if sleep_duration < 7:
                sleep_duratrion_score = 0
            elif 7 <= sleep_duration <= 9:
                sleep_duratrion_score = 5
            elif sleep_duration > 9:
                sleep_duratrion_score = 3
            sleep_list.append(duration_score_dict.get(sleep_duratrion_score))
            shift_worker_list.append(
                duration_score_dict.get(sleep_duratrion_score))

            # gastric
            # question 37 a,b,c  >1 selected heigh risk
            bloating = consumer_internal_markers.get("QUE-DET-DIG-Y-BT")
            acidic = consumer_internal_markers.get("QUE-DET-DIG-Y-AC")
            constipation = consumer_internal_markers.get("QUE-DET-DIG-Y-CN")
            if bloating == True or acidic == True or constipation == True:
                gastro_list.append("high_risk")
            else:
                gastro_list.append("low_risk")

            if "high_risk" in general_health_list or "high_risk" in weight_mgmt_list or "high_risk" in pe_list or "high_risk" in sleep_list or "high_risk" in gastro_list or "high_risk" in shift_worker_list:

                if goal_data == "General Health":
                    general_health_priority_list = [
                        "Gastro fitness", "Weight management", "Performance & endurance", "Sleep", "Shift workers"]
                    high_risk_list_count_gen_health = [gastro_list.count("high_risk"), weight_mgmt_list.count(
                        "high_risk"), pe_list.count("high_risk"), sleep_list.count("high_risk"), shift_worker_list.count("high_risk")]
                    goal = secondary_goal_priority_helper(
                        general_health_priority_list, high_risk_list_count_gen_health)

                elif goal_data == "Weight management":
                    weight_mgmt_priority_list = [
                        "General Health", "Gastro fitness", "Performance & endurance", "Sleep", "Shift workers"]
                    high_risk_list_count_weight = [general_health_list.count("high_risk"), gastro_list.count(
                        "high_risk"), pe_list.count("high_risk"), sleep_list.count("high_risk"), shift_worker_list.count("high_risk")]
                    goal = secondary_goal_priority_helper(
                        weight_mgmt_priority_list, high_risk_list_count_weight)

                elif goal_data == "Shift workers":
                    shift_worker_priority_list = [
                        "Sleep", "General Health", "Gastro fitness", "Weight management", "Performance & endurance", ]
                    high_risk_list_count_shift = [sleep_list.count("high_risk"), general_health_list.count(
                        "high_risk"), gastro_list.count("high_risk"), weight_mgmt_list.count("high_risk"), pe_list.count("high_risk")]
                    goal = secondary_goal_priority_helper(
                        shift_worker_priority_list, high_risk_list_count_shift)

                elif goal_data == "Sleep":
                    sleep_priority_list = ["Weight management", "Gastro fitness",
                                           "General Health", "Performance & endurance", "Shift workers"]
                    high_risk_list_count_sleep = [weight_mgmt_list.count("high_risk"), gastro_list.count(
                        "high_risk"), general_health_list.count("high_risk"), pe_list.count("high_risk"), shift_worker_list.count("high_risk")]
                    goal = secondary_goal_priority_helper(
                        sleep_priority_list, high_risk_list_count_sleep)

                elif goal_data == "Performance & endurance":
                    performance_priority_list = [
                        "General Health", "Sleep", "Gastro fitness", "Weight management", "Shift workers"]
                    high_risk_list_count_performance = [general_health_list.count("high_risk"), sleep_list.count(
                        "high_risk"), gastro_list.count("high_risk"), weight_mgmt_list.count("high_risk"), shift_worker_list.count("high_risk")]
                    goal = secondary_goal_priority_helper(
                        performance_priority_list, high_risk_list_count_performance)

            elif "moderate_risk" in general_health_list or "moderate_risk" in weight_mgmt_list or "moderate_risk" in pe_list or "moderate_risk" in sleep_list or "moderate_risk" in gastro_list or "moderate_risk" in shift_worker_list:

                if goal_data == "General Health":
                    general_health_priority_list = [
                        "Gastro fitness", "Weight management", "Performance & endurance", "Sleep", "Shift workers"]
                    moderate_risk_list_count_gen_health = [gastro_list.count("moderate_risk"), weight_mgmt_list.count(
                        "moderate_risk"), pe_list.count("moderate_risk"), sleep_list.count("moderate_risk"), shift_worker_list.count("moderate_risk")]
                    goal = secondary_goal_priority_helper(
                        general_health_priority_list, moderate_risk_list_count_gen_health)

                elif goal_data == "Weight management":
                    weight_mgmt_priority_list = [
                        "General Health", "Gastro fitness", "Performance & endurance", "Sleep", "Shift workers"]
                    moderate_risk_list_count_weight = [general_health_list.count("moderate_risk"), gastro_list.count(
                        "moderate_risk"), pe_list.count("moderate_risk"), sleep_list.count("moderate_risk"), shift_worker_list.count("moderate_risk")]
                    goal = secondary_goal_priority_helper(
                        weight_mgmt_priority_list, moderate_risk_list_count_weight)

                elif goal_data == "Shift workers":
                    shift_worker_priority_list = [
                        "Sleep", "General Health", "Gastro fitness", "Weight management", "Performance & endurance", ]
                    moderate_risk_list_count_shift = [sleep_list.count("moderate_risk"), general_health_list.count(
                        "moderate_risk"), gastro_list.count("moderate_risk"), weight_mgmt_list.count("moderate_risk"), pe_list.count("moderate_risk")]
                    goal = secondary_goal_priority_helper(
                        shift_worker_priority_list, moderate_risk_list_count_shift)

                elif goal_data == "Sleep":
                    sleep_priority_list = ["Weight management", "Gastro fitness",
                                           "General Health", "Performance & endurance", "Shift workers"]
                    moderate_risk_list_count_sleep = [weight_mgmt_list.count("moderate_risk"), gastro_list.count(
                        "moderate_risk"), general_health_list.count("moderate_risk"), pe_list.count("moderate_risk"), shift_worker_list.count("moderate_risk")]
                    goal = secondary_goal_priority_helper(
                        sleep_priority_list, moderate_risk_list_count_sleep)

                elif goal_data == "Performance & endurance":
                    performance_priority_list = [
                        "General Health", "Sleep", "Gastro fitness", "Weight management", "Shift workers"]
                    moderate_risk_list_count_performance = [general_health_list.count("moderate_risk"), sleep_list.count(
                        "moderate_risk"), gastro_list.count("moderate_risk"), weight_mgmt_list.count("moderate_risk"), shift_worker_list.count("moderate_risk")]
                    goal = secondary_goal_priority_helper(
                        performance_priority_list, moderate_risk_list_count_performance)

            elif "low_risk" in general_health_list or "low_risk" in weight_mgmt_list or "low_risk" in pe_list or "low_risk" in sleep_list or "low_risk" in gastro_list or "low_risk" in shift_worker_list:

                if goal_data == "General Health":
                    general_health_priority_list = [
                        "Gastro fitness", "Weight management", "Performance & endurance", "Sleep", "Shift workers"]
                    low_risk_list_count_gen_health = [gastro_list.count("low_risk"), weight_mgmt_list.count(
                        "low_risk"), pe_list.count("low_risk"), sleep_list.count("low_risk"), shift_worker_list.count("low_risk")]
                    goal = secondary_goal_priority_helper(
                        general_health_priority_list, low_risk_list_count_gen_health)

                elif goal_data == "Weight management":
                    weight_mgmt_priority_list = [
                        "General Health", "Gastro fitness", "Performance & endurance", "Sleep", "Shift workers"]
                    low_risk_list_count_weight = [general_health_list.count("low_risk"), gastro_list.count(
                        "low_risk"), pe_list.count("low_risk"), sleep_list.count("low_risk"), shift_worker_list.count("low_risk")]
                    goal = secondary_goal_priority_helper(
                        weight_mgmt_priority_list, low_risk_list_count_weight)

                elif goal_data == "Shift workers":
                    shift_worker_priority_list = [
                        "Sleep", "General Health", "Gastro fitness", "Weight management", "Performance & endurance", ]
                    low_risk_list_count_shift = [sleep_list.count("low_risk"), general_health_list.count(
                        "low_risk"), gastro_list.count("low_risk"), weight_mgmt_list.count("low_risk"), pe_list.count("low_risk")]
                    goal = secondary_goal_priority_helper(
                        shift_worker_priority_list, low_risk_list_count_shift)

                elif goal_data == "Sleep":
                    sleep_priority_list = ["Weight management", "Gastro fitness",
                                           "General Health", "Performance & endurance", "Shift workers"]
                    low_risk_list_count_sleep = [weight_mgmt_list.count("low_risk"), gastro_list.count(
                        "low_risk"), general_health_list.count("low_risk"), pe_list.count("low_risk"), shift_worker_list.count("low_risk")]
                    goal = secondary_goal_priority_helper(
                        sleep_priority_list, low_risk_list_count_sleep)

                elif goal_data == "Performance & endurance":
                    performance_priority_list = [
                        "General Health", "Sleep", "Gastro fitness", "Weight management", "Shift workers"]
                    low_risk_list_count_performance = [general_health_list.count("low_risk"), sleep_list.count(
                        "low_risk"), gastro_list.count("low_risk"), weight_mgmt_list.count("low_risk"), shift_worker_list.count("low_risk")]
                    goal = secondary_goal_priority_helper(
                        performance_priority_list, low_risk_list_count_performance)
            return goal

    @staticmethod
    def get_Contributors_fiber_intake(npe_answers):
        if round(npe_answers["QUE-HEH-GR-PER"]) and round(npe_answers["QUE-HEH-BL-PER"]) and round(npe_answers["QUE-HEH-VF-PER"]):
            intake="The contributors of your fiber intake are cereals, Legumes, vegetables and fruits only."
        elif round(npe_answers["QUE-HEH-GR-PER"]) and round(npe_answers["QUE-HEH-BL-PER"]):
            intake="The contributors of your fiber intake are cereals and Legumes only."
        elif round(npe_answers["QUE-HEH-BL-PER"]) and round(npe_answers["QUE-HEH-VF-PER"]):
            intake="The contributors of your fiber intake are Legumes, vegetables and fruits only."
        elif round(npe_answers["QUE-HEH-GR-PER"]) and round(npe_answers["QUE-HEH-VF-PER"]):
            intake="The contributors of your fiber intake are cereals, vegetables and fruits only."
        elif round(npe_answers["QUE-HEH-GR-PER"]):
            intake="The contributors of your fiber intake is cereals only."
        elif round(npe_answers["QUE-HEH-BL-PER"]):
            intake="The contributors of your fiber intake is Legumes only."
        elif round(npe_answers["QUE-HEH-VF-PER"]):
            intake="The contributors of your fiber intake are vegetables and fruits only."
        else:
            intake="The contributors of your fiber intake is None."
        return intake 

