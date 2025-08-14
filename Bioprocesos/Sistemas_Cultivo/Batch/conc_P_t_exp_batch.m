function [t_valores, P_valores] = conc_P_t_exp_batch(u_max, y_xs, y_ps, P_inicial, X_inicial, t_final, t_lag, delta_t)
  % Modelo de crecimiento exponencial para la biomasa. [P](t) = [X](0) * (Y_p/s)/ (Y_x/s)*(e^(u_max*(t-t_lag)-1) + [P](0) en sistema Batch.
  % En Batch se cumple d[P]/dt = rp = qp*[X](t) =(u_max*y_ps/y_xs)*[X](0)*e^(u_max*(t-t_lag)
  % Parámetros: u_max (velocidad de crecimiento específica máxima), y_xs (rendimiento celular),y_ps (rendimiento del producto), P_inicial (producto inicial), X_inicial (inoculo), t_final (tiempo final estimado), t_lag (tiempo de la fase lag), delta_t (separación entre puntos de tiempo).
  % Salida : Matriz rectangular de valores de tiempo y valores de producto.

  if nargin < 7, t_lag = 0; endif %Si no se especifíca, la fase lag es nula.
  if nargin < 8, delta_t = 0.1; endif

  t_valores = [];
  P_valores = [];
  t_valores_lag = [];
  P_valores_lag = [];

  if t_lag > 0
    t_valores_lag = 0: delta_t : t_lag;
    N = numel(t_valores_lag);
    P_valores_lag = P_inicial*ones(1,N);
  endif

  t_valores_exp = t_lag : delta_t : t_final;
  N = numel(t_valores_exp);
  P_valores_exp = zeros(1, N);
  P_valores_exp(1) = P_inicial;
  for k = 2:N
    P_valores_exp(k) = P_inicial + (X_inicial*y_ps/y_xs)*(exp(u_max*(t_valores_exp(k)-t_lag))-1);
  endfor

  P_valores = [P_valores_lag(:); P_valores_exp(:)];
  t_valores = [t_valores_lag(:); t_valores_exp(:)];

endfunction
