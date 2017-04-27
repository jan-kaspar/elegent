FROM base/archlinux

RUN pacman --noconfirm -Sy gcc make cmake git gsl libglvnd openssl

RUN for p in \
		hepmc-2.06.09-1-x86_64.pkg.tar.xz \
		libafterimage-1.20-2-x86_64.pkg.tar.xz \
		unuran-1.8.1-1-x86_64.pkg.tar.xz \
		root-6.08.06-1-x86_64.pkg.tar.xz; \
	do \
		curl "http://jkaspar.web.cern.ch/jkaspar/share/archlinux/$p" > "$p"; \
		pacman --noconfirm -U "$p"; \
	done
