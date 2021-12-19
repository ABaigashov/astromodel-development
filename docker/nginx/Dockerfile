FROM nginx

COPY ${PWD}/docker/nginx/nginx.conf /etc/nginx/conf.d/default.conf
COPY ${PWD}/server/www /www

COPY ${PWD}/configurator /configurator

ENTRYPOINT ["nginx", "-g", "daemon off;"]