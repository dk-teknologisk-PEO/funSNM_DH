function [long, lat] = city_position(name)
    data = webread('https://api.dataforsyningen.dk/steder?hovedtype=Bebyggelse&undertype=by&','primærtnavn',name);
    position = data.visueltcenter;
    long=position(1);
    lat = position(2);
end