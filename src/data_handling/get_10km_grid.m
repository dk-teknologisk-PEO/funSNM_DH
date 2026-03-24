function grid_id_10km = get_10km_grid(long, lat)
    data = webread('https://api.dataforsyningen.dk/adgangsadresser/reverse?','x',long,'y',lat);
    grid = data.DDKN;
    grid_id_10km = grid.km10;
    disp(grid_id_10km)
end