#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 20:12:01 2020

@author: yangwenli
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import utils

class CTAPortfolio:
    
    def __init__(self, hist_data, loss_stop = -0.05, profit_stop = None, level = 1, transaction_cost = 0.0005, long_open_signal = None, long_close_signal = None, short_open_signal = None, short_close_signal = None):
        self._hist_data = hist_data
        self._returns = None
        self._loss_stop = loss_stop
        self._profit_stop = profit_stop
        self._returns_list = None
        self._level = level
        self._transaction_cost = transaction_cost
        
        if long_open_signal is not None:
            if long_close_signal is not None:
                self._buy_open_times = pd.Series(long_open_signal.loc[long_open_signal].index)
                self._buy_close_times = pd.Series(long_close_signal.loc[long_close_signal].index)
            else:
                raise ValueError('买入开仓函数必须有平仓函数相匹配！')
        else:
            self._buy_open_times = self._buy_close_times = None
            
        if short_open_signal is not None:
            if short_close_signal is not None:
                self._sell_open_times = pd.Series(short_open_signal.loc[short_open_signal].index)
                self._sell_close_times = pd.Series(short_close_signal.loc[short_close_signal].index)
            else:
                raise ValueError('卖出开仓函数必须有平仓函数相匹配！')
        else:
            self._sell_open_times = self._sell_close_times = None
                
    def _cal_returns(self):
        if self._returns_list is None:
            self._find_open_period(self._buy_open_times, self._buy_close_times, self._sell_open_times, self._sell_close_times)
        
        returns = []
        for close_series in self._returns_list['long']:
            period_returns = close_series.pct_change().dropna()
            #扣除手续费
            period_returns.iloc[0] = (period_returns.iloc[0] + 1) / (1 + self._transaction_cost) - 1
            period_returns.iloc[-1] = (period_returns.iloc[-1] + 1) * (1 - self._transaction_cost) - 1
            returns.append(period_returns)
        
        for close_series in self._returns_list['short']:
            period_returns = close_series.pct_change().dropna()
            #扣除手续费
            period_returns.iloc[0] = (period_returns.iloc[0] + 1) / (1 + self._transaction_cost) - 1
            period_returns.iloc[-1] = (period_returns.iloc[-1] + 1) * (1 - self._transaction_cost) - 1
            returns.append(period_returns)
        returns = pd.concat(returns)
        returns = returns.reindex(self._hist_data.index).fillna(0)
        return returns
    
    @property
    def returns(self):
        return self._cal_returns()
    
    @property
    def long_returns(self):
        if self._returns_list is None:
            self._find_open_period(self._buy_open_times, self._buy_close_times, self._sell_open_times, self._sell_close_times)
        
        returns = []
        for close_series in self._returns_list['long']:
            period_returns = close_series.pct_change().dropna()
            #扣除手续费
            period_returns.iloc[0] = (period_returns.iloc[0] + 1) / (1 + self._transaction_cost) - 1
            period_returns.iloc[-1] = (period_returns.iloc[-1] + 1) * (1 - self._transaction_cost) - 1
            returns.append(period_returns)
            
        returns = pd.concat(returns)
        returns = returns.reindex(self._hist_data.index).fillna(0)
        return returns
    
    @property
    def short_returns(self):
        if self._returns_list is None:
            self._find_open_period(self._buy_open_times, self._buy_close_times, self._sell_open_times, self._sell_close_times)
        
        returns = []
        for close_series in self._returns_list['short']:
            period_returns = close_series.pct_change().dropna()
            #扣除手续费
            period_returns.iloc[0] = (period_returns.iloc[0] + 1) / (1 + self._transaction_cost) - 1
            period_returns.iloc[-1] = (period_returns.iloc[-1] + 1) * (1 - self._transaction_cost) - 1
            returns.append(period_returns)
            
        returns = pd.concat(returns)
        returns = returns.reindex(self._hist_data.index).fillna(0)
        return returns
    
    def _find_open_period(self, long_open_times = None, long_close_times = None, short_open_times = None, short_close_times = None, direction = 'long'):
        '''找到开仓和平仓之间的时间段，只计算相应时间段内收益率即可'''
        times_list = [times if times is not None else pd.Series() for times in [long_open_times, long_close_times, short_open_times, short_close_times]]
        returns_list = {'long':[], 'short':[]}
        while times_list[0].shape[0] != 0 and times_list[2].shape[0] != 0:
            start_time = min(times_list[0].iloc[0], times_list[2].iloc[0])
            if start_time == times_list[0].iloc[0]:
                #如果开了多单
                direction = 1
            else:
                #如果开了空单
                direction = -1
            
            times_list, period_close = self._find_close_series(times_list, start_time, direction)
            
            if direction == 1:
                returns_list['long'].append(period_close)
            else:
                returns_list['short'].append(period_close)
            
        
        while times_list[0].shape[0] != 0:
            start_time = times_list[0].iloc[0]
            direction = 1
            times_list, period_close = self._find_close_series(times_list, start_time, direction)
            returns_list['long'].append(period_close)
            
        while times_list[2].shape[0] != 0:
            start_time = times_list[2].iloc[0]
            direction = -1
            times_list, period_close = self._find_close_series(times_list, start_time, direction)
            returns_list['short'].append(period_close)
        
        self._returns_list = returns_list
            
    
    def _find_close_series(self, times_list, start_time, direction):
        if direction == 1:
            #开多单
            times_list[1] = times_list[1][times_list[1] > start_time]
            if len(times_list[1]) != 0:
                end_time = times_list[1].iloc[0]
            else:
                end_time = self._hist_data.index[-1]
        else:
            times_list[3] = times_list[3][times_list[3] > start_time]
            if len(times_list[3]) != 0:
                end_time = times_list[3].iloc[0]
            else:
                end_time = self._hist_data.index[-1]
                
        period_close = self._hist_data['close'][start_time:end_time]
        period_close = period_close.iloc[0] + (period_close - period_close.iloc[0]) * direction
        
        if self._loss_stop is not None:
            below_loss_stop = period_close / period_close.iloc[0] - 1 <= self._loss_stop
            if below_loss_stop.any():
                below_loss_time = below_loss_stop.idxmax()
                period_close = period_close[:below_loss_time]
            # max_value = period_close.iloc[0]
            # for i in range(len(period_close)):
            #     value_now = period_close.iloc[i]
            #     max_value = max(max_value, value_now)
            #     if value_now / max_value - 1 <= self._loss_stop:
            #         break
            # period_close = period_close.iloc[:(i + 1)]
                    
        if self._profit_stop is not None:
            above_profit_stop = period_close / period_close.iloc[0] - 1 >= self._profit_stop
            if above_profit_stop.any():
                above_profit_time = above_profit_stop.idxmax()
                period_close = period_close[:above_profit_time]
            
        stop_time = period_close.index[-1]
        times_list = [times[times >= stop_time] if len(times) > 0 else times for times in times_list]
        
        return times_list, period_close
    
    def _function_to_signals(self, funcs):
        '''传入函数的字典计算信号'''
        if not isinstance(funcs, tuple):
            funcs = ('and', funcs)
            
        combine_func = utils.and_ if funcs[0] == 'and' else utils.or_
        result_list = []
        for i in range(1, len(funcs)):
            resulti = self._function_to_signals(funcs[i]) if isinstance(funcs[i], tuple) else funcs[i](self._hist_data)
            result_list.append(resulti)
        signals = combine_func(*result_list)
        return signals
    
    def _signals_time(self, signals):    
        signals_time = pd.Series(signals.index[signals])
        return signals_time
    
    def buy_open_signals(self, buy_open_funcs):
        '''计算买入开仓信号'''
        return self._signals_time(self._function_to_signals(buy_open_funcs))
    
    def buy_close_signals(self, buy_close_funcs):
        '''计算买入平仓信号'''
        return self._signals_time(self._function_to_signals(buy_close_funcs))
    
    def sell_open_signals(self, sell_open_funcs):
        '''计算买入开仓信号'''
        return self._signals_time(self._function_to_signals(sell_open_funcs))
    
    def sell_close_signals(self, sell_close_funcs):
        '''计算买入平仓信号'''
        return self._signals_time(self._function_to_signals(sell_close_funcs))
    
    @property
    def trade_count(self):
        '''按年统计交易次数'''
        long_open_time = dict((returns.index[0], 1) for returns in self._returns_list['long'])
        short_open_time = dict((returns.index[0], 1) for returns in self._returns_list['short'])
        long_open_time = pd.Series(long_open_time)
        short_open_time = pd.Series(short_open_time)
        combine_open_time = pd.concat([long_open_time, short_open_time], axis = 1)
        combine_open_time.columns = ['long', 'short']
        combine_open_time.index = pd.to_datetime(combine_open_time.index)
        return combine_open_time.resample('Y').sum(), combine_open_time.sum()
    
    @property
    def holding_days(self):
        '''统计每年持有自然日'''
        long_holding_days = dict((returns.index[0], (returns.index[-1] - returns.index[0]).days + 1) for returns in self._returns_list['long'])
        short_holding_days = dict((returns.index[0], (returns.index[-1] - returns.index[0]).days + 1) for returns in self._returns_list['short'])
        long_holding_days = pd.Series(long_holding_days)
        short_holding_days = pd.Series(short_holding_days)
        combine_holding_days = pd.concat([long_holding_days, short_holding_days], axis = 1)
        combine_holding_days.columns = ['long', 'short']
        combine_holding_days.index = pd.to_datetime(combine_holding_days.index)
        return combine_holding_days.resample('Y').sum(), combine_holding_days.sum()
    
    @property
    def win_ratio(self):
        '''统计胜率'''
        # long_profit_time = dict((close.index[0], close.iloc[-1] > close.iloc[0]) for close in self._returns_list['long'])
        # short_profit_time = dict((close.index[0], close.iloc[-1] > close.iloc[0]) for close in self._returns_list['short'])
        # long_profit_time = pd.Series(long_profit_time)
        # short_profit_time = pd.Series(short_profit_time)
        # long_win_ratio = long_profit_time.mean()
        # short_win_ratio = short_profit_time.mean()
        # long_win_yearly_ratio = long_profit_time.resample('Y').mean()
        # short_win_yearly_ratio = short_profit_time.resample('Y').mean()
        # win_ratio = pd.Series({'long':long_win_ratio, 'short':short_win_ratio})
        # yearly_win_ratio = pd.concat([long_win_yearly_ratio, short_win_yearly_ratio], axis = 1)
        # yearly_win_ratio.columns = ['long', 'short']
        # return win_ratio, yearly_win_ratio
        long_win_time = dict((returns.index[0], returns.iloc[-1] > returns.iloc[0]) for returns in self._returns_list['long'])
        short_win_time = dict((returns.index[0], returns.iloc[-1] > returns.iloc[0]) for returns in self._returns_list['short'])
        long_win_time = pd.Series(long_win_time)
        short_win_time = pd.Series(short_win_time)
        combine_win_time = pd.concat([long_win_time, short_win_time], axis = 1).replace(np.nan, False)
        combine_win_time.columns = ['long', 'short']
        combine_win_time.index = pd.to_datetime(combine_win_time.index)
        return combine_win_time.resample('Y').mean(), combine_win_time.mean()
    
    @property
    def win_count(self):
        '''统计盈利次数'''
        long_win_time = dict((returns.index[0], returns.iloc[-1] > returns.iloc[0]) for returns in self._returns_list['long'])
        short_win_time = dict((returns.index[0], returns.iloc[-1] > returns.iloc[0]) for returns in self._returns_list['short'])
        long_win_time = pd.Series(long_win_time)
        short_win_time = pd.Series(short_win_time)
        combine_win_time = pd.concat([long_win_time, short_win_time], axis = 1)
        combine_win_time.columns = ['long', 'short']
        combine_win_time.index = pd.to_datetime(combine_win_time.index)
        return combine_win_time.resample('Y').sum(), combine_win_time.sum()
    
    @property
    def annual_returns(self):
        returns = self.returns.resample('D').apply(lambda s:(1 + s).cumprod().iloc[-1] - 1 if len(s) > 0 else np.nan).dropna()
        return returns.resample('Y').mean() * 244, returns.mean() * 244
    
    @property
    def profit_loss_ratio(self):
        '''统计盈亏比'''
        long_open_returns = dict((close.index[0], close.iloc[-1] / close.iloc[0] - 1) for close in self._returns_list['long'])
        short_open_returns = dict((close.index[0], close.iloc[-1] / close.iloc[0] - 1) for close in self._returns_list['short'])
        long_open_returns = pd.Series(long_open_returns)
        short_open_returns = pd.Series(short_open_returns)
        
        long_profit_loss_ratio = -long_open_returns[long_open_returns > 0].mean() / long_open_returns[long_open_returns < 0].mean()
        short_profit_loss_ratio = -short_open_returns[short_open_returns > 0].mean() / short_open_returns[short_open_returns < 0].mean()     
        
        open_returns = pd.concat([long_open_returns, short_open_returns])
        profit_loss_ratio = -open_returns[open_returns > 0].mean() / open_returns[open_returns < 0].mean()
        
        if len(long_open_returns) != 0:
            yearly_long_profit_loss_ratio = -long_open_returns[long_open_returns > 0].resample('Y').mean().fillna(0) / long_open_returns[long_open_returns < 0].resample('Y').mean().fillna(0)
        else:
            yearly_long_profit_loss_ratio = long_open_returns
        
        if len(short_open_returns) != 0:
            yearly_short_profit_loss_ratio = -short_open_returns[short_open_returns > 0].resample('Y').mean().fillna(0) / short_open_returns[short_open_returns < 0].resample('Y').mean().fillna(0)
        else:
            yearly_short_profit_loss_ratio = short_open_returns
        
        yearly_profit_loss_ratio = pd.concat([yearly_long_profit_loss_ratio, yearly_short_profit_loss_ratio], axis = 1)
        yearly_profit_loss_ratio.columns = ['long', 'short']
        
        return yearly_profit_loss_ratio, pd.Series({'long':long_profit_loss_ratio,'short':short_profit_loss_ratio}), profit_loss_ratio
    
    @property
    def annual_vol(self):
        returns = self.returns.resample('D').apply(lambda s:(1 + s).cumprod().iloc[-1] - 1 if len(s) > 0 else np.nan).dropna()
        return returns.resample('Y').std() * 244 ** 0.5, returns.std() * 244 ** 0.5
    
    def _max_draw_down(self, s):
        '''计算最大回撤的函数'''
        s = (s + 1).cumprod()
        max_value = 1
        mdd = 0
        for value in s:
            draw_down = value / max_value - 1
            if draw_down < mdd: mdd = draw_down
            if value > max_value: max_value = value
        return mdd
    
    @property
    def period_net_value(self):
        returns = self.returns
        return returns.resample('Y').apply(lambda s:(1 + s).cumprod().iloc[-1]), (returns + 1).cumprod().iloc[-1]
    
    @property
    def max_drawdown(self):
        returns = self.returns
        return returns.resample('Y').apply(self._max_draw_down), self._max_draw_down(returns)
    
    @property
    def sharpe_ratio(self):
        return self.annual_returns[0] / self.annual_vol[0], self.annual_returns[1] / self.annual_vol[1]
    
    @property
    def information_ratio(self):
        returns = pd.concat([self.returns, self._hist_data['returns']], axis = 1).fillna(0)
        returns.columns = ['portfolio', 'benchmark']
        excess_returns = returns['portfolio'] - returns['benchmark']
        return excess_returns.resample('Y').mean() / excess_returns.resample('Y').std() * 244 ** 0.5, excess_returns.mean() / excess_returns.std() * 244 ** 0.5
    
    def plot_net_value(self, subset = 'long', title = None, start_date = None):
        if subset == 'long':
            returns = self.long_returns
        elif subset == 'short':
            returns = self.short_returns
        else:
            returns = self.returns
        
        returns = pd.concat([returns, self._hist_data.returns], axis = 1).fillna(0)
        returns.columns = ['portfolio', 'benchmark']
        if start_date is None:
            return (1 + returns).cumprod().plot(title = title)
        else:
            return (1 + returns[start_date:]).cumprod().plot(title = title)
    
    def summary(self, yearly = False, title = None, savefig = False, plot = True):
        if self._returns_list is None:
            self._find_open_period(self._buy_open_times, self._buy_close_times, self._sell_open_times, self._sell_close_times)
        if plot:
            self.plot_net_value(subset = 'long-and-short', title = title)
        if savefig:
            plt.savefig(f'{title}.png')
        plt.close()
        
        # summary_dict = {}
        # summary_dict['年化收益率'] = '{:.2%}'.format(self.annual_returns)
        # summary_dict['年化波动率'] = '{:.2%}'.format(self.annual_vol)
        # summary_dict['夏普比率'] = '{:.2}'.format(self.sharpe_ratio)
        # summary_dict['盈亏比'] = '{:.2}'.format(self.profit_loss_ratio[1])
        # summary_dict['最大回撤'] = '{:.2%}'.format(self.max_drawdown)
        # summary_dict['开仓次数'] = int(self.trade_count()[1].sum())
        
        # return pd.Series(summary_dict)
        
        statistics = [self.trade_count, self.holding_days, self.win_count, self.win_ratio, self.profit_loss_ratio,
                      self.period_net_value, self.annual_returns, self.annual_vol, self.sharpe_ratio, self.information_ratio, self.max_drawdown]
        columns = ['开仓次数', '持有交易日', '盈利次数', '胜率', '盈亏比', '区间净值', '年化收益率', '年化波动率', '收益风险比', '信息比率', '最大回撤']
        
        columns2 = []
        for i in range(len(columns)):
            c = columns[i]
            if i <= 4:
                columns2 += [f'{c}多头', f'{c}空头']
            else:
                columns2 += [c]
        if yearly:
            summaries = pd.concat([s[0] for s in statistics], axis = 1)
            summaries.columns = columns2
        else:
            summaries = pd.Series([s[1] if isinstance(s[1], pd.Series) else s[0] for s in statistics], index = columns2)
        
        return summaries