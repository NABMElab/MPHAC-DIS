function encoded = rsEncoder(msg, m, prim_poly, n, k)
    % RSENCODER Encode message with the Reed-Solomon algorithm
    % m is the number of bits per symbol
    % prim_poly: Primitive polynomial p(x). Ie for DM is 301
    % k is the size of the message
    % n is the total size (k+redundant)
    % Example: msg = uint8('Test')
    % enc_msg = rsEncoder(msg, 8, 301, 12, numel(msg));

    % Get the alpha
    alpha = gf(2, m, prim_poly);

    % Get the Reed-Solomon generating polynomial g(x)
    g_x = genpoly(k, n, alpha);

    % Multiply the information by X^(n-k), or just pad with zeros at the end to
    % get space to add the redundant information
    msg_padded = gf([msg zeros(1, n - k)], m, prim_poly);

    % Get the remainder of the division of the extended message by the
    % Reed-Solomon generating polynomial g(x)
    [~, remainder] = deconv(msg_padded, g_x);

    % Now return the message with the redundant information
    encoded = msg_padded - remainder;

end

% Find the Reed-Solomon generating polynomial g(x), by the way this is the
% same as the rsgenpoly function on matlab
function g = genpoly(k, n, alpha)
    g = 1;
    % A multiplication on the galois field is just a convolution
    for k = mod(1 : n - k, n)
        g = conv(g, [1 alpha .^ (k)]);
    end
end