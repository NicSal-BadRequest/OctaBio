function [t_valores, X_valores] = conc_X_t_exp_batch(u_max, X_inicial, t_final, t_lag, delta_t)

  % Modelo de crecimiento exponencial para la biomasa. [X](t) = [X]0 e^(u_max*(t-t_lag)). Considera la fase lag y sistema Batch ([S] >> Ks de modo que u = u_max).
  % Parámetros: u_max (velocidad de crecimiento específica máxima), x_inicial (biomasa inicial o inoculo), t_final (tiempo final estimado), t_lag (tiempo de la fase lag), delta_t (separación entre puntos de tiempo).
  % Salida : Matriz rectangular de valores de tiempo y valores de biomasa.

  if nargin < 4, t_lag = 0; endif %Si no se especifíca, la fase lag es nula.
  if nargin < 5, delta_t = 0.1; endif

  t_valores = [];
  X_valores = [];
  t_valores_lag = [];
  X_valores_lag = [];

  if t_lag > 0
    t_valores_lag = 0:delta_t:t_lag;
    N = numel(t_valores_lag);
    X_valores_lag = X_inicial*ones(1,N);
  endif

   t_valores_exp = t_lag : delta_t : t_final;
   N = numel(t_valores_exp);
   X_valores_exp = zeros(1, N);
   X_valores_exp(1) = X_inicial;
   for k = 2:N
    X_valores_exp(k) = X_inicial*exp(u_max*(t_valores_exp(k)-t_lag));
   endfor


  X_valores = [X_valores_lag(:); X_valores_exp(:)];
  t_valores = [t_valores_lag(:); t_valores_exp(:)];

endfunction

