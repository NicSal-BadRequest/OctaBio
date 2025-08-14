function [t_valores, S_valores] = conc_S_t_exp_batch(u_max, y_xs, S_inicial, X_inicial, t_final, t_lag, delta_t)
  % Modelo de crecimiento exponencial para la biomasa. [S](t) = [X](0)/ (Y_x/s)*(1-e^(u_max*(t-t_lag)) + [S] (0) en sistema Batch.
  % En Batch se cumple d[S]/dt = -rs = -qs*[X](t) =-(u_max/y_xs)*[X](0)*e^(u_max*(t-t_lag)
  % Parámetros: u_max (velocidad de crecimiento específica máxima), y_xs (rendimiento celular), S_inicial (sustrato inicial), X_inicial (inoculo), t_final (tiempo final estimado), t_lag (tiempo de la fase lag), delta_t (separación entre puntos de tiempo).
  % Salida : Matriz rectangular de valores de tiempo y valores de sustrato.

  if nargin < 6, t_lag = 0; endif %Si no se especifíca, la fase lag es nula.
  if nargin < 7, delta_t = 0.1; endif

  t_valores = [];
  S_valores = [];
  t_valores_lag = [];
  S_valores_lag = [];

  if t_lag > 0
    t_valores_lag = 0: delta_t : t_lag;
    N = numel(t_valores_lag);
    S_valores_lag = S_inicial*ones(1,N);
  endif

  t_valores_exp = t_lag : delta_t : t_final;
  N = numel(t_valores_exp);
  S_valores_exp = zeros(1, N);
  S_valores_exp(1) = S_inicial;
  for k = 2:N
    S_valores_exp(k) = S_inicial + (X_inicial/y_xs)*(1-exp(u_max*(t_valores_exp(k)-t_lag)));
    if (S_valores_exp(k) < 0) %No puede haber sustrato negativo.
      t_valores_exp = t_valores_exp(1:k-1);
      S_valores_exp = S_valores_exp(1:k-1);
      break;
    endif
  endfor


  S_valores = [S_valores_lag(:); S_valores_exp(:)];
  t_valores = [t_valores_lag(:); t_valores_exp(:)];

endfunction
