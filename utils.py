#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 22:44:00 2020

@author: yangwenli
"""

import pandas as pd

def read_data(file_path):
    data = pd.read_csv(file_path)
    data['date'] = pd.to_datetime(data['date'].apply(str))
    if 'code' in data.columns:
        data = data.set_index(['code', 'date'])
    else:
        data = data.set_index('date')
    return data

def and_(*conds):
    '''计算pandas.Series格式的条件的和'''
    res = conds[0]
    for cond in conds[1:]:
        res = res & cond
    return res

def or_(*conds):
    '''计算pandas.Series格式的条件的或'''
    res = conds[0]
    for cond in conds[1:]:
        res = res | cond
    return res

def get_dominate_data(hist_data, price_columns = ['open', 'high', 'close', 'low']):
    '''获取期货合约的主力合约数据'''
    #每个时刻成交量最大的合约作为下一时刻的主力合约
    dominate_data_index = hist_data['volume'].groupby('date').idxmax()
    dominate_data_index = dominate_data_index.apply(lambda t:t[0])
    dominate_data_index.name = 'dominate_code'
    last_dominate_data_index = dominate_data_index.shift(1)
    last_dominate_data_index.name = 'last_dominate_code'
    next_dominate_data_index = dominate_data_index.shift(1)
    next_dominate_data_index.name = 'next_dominate_code'
    dominate_data_index = dominate_data_index.reset_index().set_index(['dominate_code', 'date']).index
    last_dominate_data_index = last_dominate_data_index.dropna().reset_index().set_index(['last_dominate_code', 'date']).index    
    next_dominate_data_index = next_dominate_data_index.dropna().reset_index().set_index(['next_dominate_code', 'date']).index
    #未考虑合约一上市就变为主力合约或一直为主力合约直到退市的情况
    sub_data = hist_data[price_columns]
    dominate_sub_data = sub_data.reindex(dominate_data_index).droplevel(level = 'dominate_code')
    last_dominate_data = sub_data.reindex(last_dominate_data_index).droplevel(level = 'last_dominate_code')
    next_dominate_data = sub_data.reindex(next_dominate_data_index).droplevel(level = 'next_dominate_code')
    adj_factor = last_dominate_data['close'] / dominate_sub_data['close']
    adj_factor = adj_factor.cumprod()
    for c in price_columns:
        dominate_sub_data[c] = adj_factor.shift(1) * next_dominate_data[c]
    
    codes = [x[0] for x in dominate_data_index]
    dates = [x[1] for x in dominate_data_index]
    dominate_data_index = pd.DataFrame({'code':codes, 'date':dates}).set_index('date')['code'].shift(1)
    return dominate_sub_data, dominate_data_index